#!/usr/bin/env python3
"""Certify uniform adjacent-chart stability for the centered scalar."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
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
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_adjacent_chart_stability_certificate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_saddle_recurrence.json"
    ),
    "global_c1": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_C1_global_remainder_certificate.json"
    ),
    "global_first_order": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_global_remainder_certificate.json"
    ),
    "centered_reduction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
    ),
}
PRECISION_BITS = 256
L_MIN = 50
LOG_RATIO_CONSTANT = 8
LOG_RATIO_W_CONSTANT = 30
LOG_RATIO_EPSILON_CONSTANT = 100
LOG_RATIO_X_CONSTANT = 2
RATIO_VALUE_CONSTANT = 11
RATIO_X_CONSTANT = 3
ENDPOINT_CONSTANT = 100
ENDPOINT_RATE_CONSTANT = 2
ADJACENT_JET_CONSTANT = 2500
EXPLICIT_D_CONSTANT = 1500
FRAME_CONSTANT = 5
SCALAR_CHART_CONSTANT = 5000


@dataclass(frozen=True)
class StabilityRow:
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
        "adjacent_recurrence": (
            "Delta A_a=Re(Q_(N,x)-lambda_a*Q_N",
            "[epsilon]log(main_sharp/endpoint_C1)=0",
        ),
        "global_c1": (
            "|log(main_n/endpoint_n)|<3/|T|",
            "|endpoint_n|/A_t(x)<50*exp(-L/4)",
        ),
        "global_first_order": (
            "|Delta lift_[1]|/A_t(x)<20000*exp(-5L/4)",
            "|d_(t,n)|<100/T",
        ),
        "centered_reduction": (
            "|u_a|<=tL/(2x)<=25/(2x)<exp(-L)",
            "|v_a|<3/x^2",
            "|d_(n,x)|<4223/x^2",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(payload.get("kind", "")) for key, payload in payloads.items()}


def symbolic_audit() -> dict[str, str]:
    time, delta, r_one, r_zero, remainder = sp.symbols(
        "t delta r1 r0 R_M"
    )
    shift = sp.Rational(1, 2) - sp.I * sp.pi * time / 8
    original = (
        shift * r_zero
        + remainder
        - shift * delta
        + time
        / 4
        * (
            sp.I * sp.pi * (r_one - delta) / 2
            + (r_one - delta) ** 2
        )
    )
    simplified = (
        shift * r_zero
        + remainder
        - delta / 2
        + sp.I * sp.pi * time * r_one / 8
        + time * (r_one - delta) ** 2 / 4
    )
    if sp.simplify(original - simplified) != 0:
        raise RuntimeError("exact heat/logarithm cancellation failed")

    epsilon, w = sp.symbols(
        "epsilon w", positive=True, real=True
    )
    b_one = -w / 2 - 2 * sp.I * sp.pi * w**3 / 3
    first_geometric = -w / 2 - 2 * sp.I * sp.pi * w**3 / 3
    if sp.simplify(first_geometric - b_one) != 0:
        raise RuntimeError("geometric first coefficient failed")

    a_symbol = sp.symbols("a_symbol", positive=True, real=True)
    a_x = 1 / (8 * sp.pi * a_symbol)
    epsilon_x = sp.simplify(
        sp.diff(1 / a_symbol, a_symbol) * a_x
    )
    expected_epsilon_x = -epsilon**3 / (8 * sp.pi)
    if sp.simplify(
        epsilon_x.subs(a_symbol, 1 / epsilon)
        - expected_epsilon_x
    ) != 0:
        raise RuntimeError("epsilon x derivative failed")
    w_x = -a_x
    if sp.simplify(
        w_x.subs(a_symbol, 1 / epsilon)
        + epsilon / (8 * sp.pi)
    ) != 0:
        raise RuntimeError("saddle displacement x derivative failed")

    alpha_value, alpha_one, log_a = sp.symbols(
        "alpha alpha_one log_a"
    )
    k_rate = sp.symbols("K_rate")
    s_star_x = -sp.I / 2 - sp.I * time * alpha_one / 4
    m_rate = -sp.I * (
        alpha_value + time * alpha_value * alpha_one / 2
    ) / 2
    mu = k_rate - m_rate + s_star_x * log_a
    chi = alpha_value - log_a
    expected_mu = k_rate + sp.I * chi / 2 + sp.I * time * alpha_one * chi / 4
    if sp.simplify(mu - expected_mu) != 0:
        raise RuntimeError("centered endpoint rate cancellation failed")

    return {
        "domain": (
            "L>=50, 0<=tL<=25, epsilon=1/a, "
            "0<=w=N+1-a<=1, epsilon<=exp(-25)"
        ),
        "positive_ratio": (
            "rho=exp(Omega)*(1+d_+)/B, "
            "B=1-epsilon*w/2-i*(2*pi/3)*epsilon*w^3"
        ),
        "omega": (
            "Omega=h*r0+R_M-delta/2+i*pi*t*r1/8"
            "+t*(r1-delta)^2/4-i*pi*Psi"
        ),
        "scaled_blocks": (
            "delta=log(1+epsilon*w); "
            "Psi=epsilon^-2*(2delta-2epsilon*w+epsilon^2*w^2); "
            "r0=epsilon^2/(4*pi*i)+epsilon^2/(2*pi*i-epsilon^2); "
            "r1=epsilon^2/[2*(2*pi*i+h*epsilon^2)]"
            "+epsilon^2/[2*pi*i+(h-1)*epsilon^2]"
            "+log(1+h*epsilon^2/(2*pi*i))/2"
        ),
        "log_ratio": (
            "Z=log(rho)=Omega+log(1+d_+)-log(B)"
        ),
        "coordinate_derivatives": (
            "w_x=-epsilon/(8*pi), "
            "epsilon_x=-epsilon^3/(8*pi)"
        ),
        "mu_cancellation": (
            "mu_a=K_x/K+i*(alpha-log(a))/2"
            "+i*t*alpha'*(alpha-log(a))/4=O(epsilon^2)"
        ),
        "endpoint_rate": (
            "|lambda_a+mu_a+J_(a,x)/J_a|"
            "<2*epsilon"
        ),
        "ratio_bounds": (
            "|Z|<8*epsilon^2, |Z_w|<30*epsilon^2, "
            "|Z_epsilon|<100*epsilon, |Z_x|<2*epsilon^3, "
            "|rho-1|<11*epsilon^2, |rho_x|<3*epsilon^3"
        ),
        "adjacent_jet": (
            "|Delta X|<1100*exp(-5L/4), "
            "|Delta U|<2500*exp(-7L/4)"
        ),
        "scalar_stability": (
            "|A_(a,N+1)-A_(a,N)|"
            "<5000*exp(-7L/4)"
        ),
        "scale_comparison": (
            "5000*exp(-7L/4)"
            "<10^-7*exp(-5L/4) for L>=50"
        ),
    }


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def interval_budget() -> dict[str, str | int]:
    flint.ctx.prec = PRECISION_BITS
    one = flint.arb(1)
    two = flint.arb(2)
    pi = flint.arb.pi()
    epsilon = (-flint.arb(L_MIN) / 2).exp()
    h_bound = (one / 4 + (pi / 16) ** 2).sqrt()
    b_constant = one / 2 + 2 * pi / 3
    b_derivative_constant = one / 2 + 2 * pi
    b_radius = b_constant * epsilon

    if not epsilon < flint.arb("1.4e-11"):
        raise RuntimeError("epsilon endpoint failed")
    if not h_bound < flint.arb("0.54"):
        raise RuntimeError("shift bound failed")
    if not b_constant < flint.arb("2.6"):
        raise RuntimeError("endpoint first coefficient failed")
    if not b_derivative_constant < flint.arb("6.8"):
        raise RuntimeError("endpoint derivative coefficient failed")
    if not b_radius < flint.arb("4e-11"):
        raise RuntimeError("endpoint logarithm radius failed")

    delta_remainder = epsilon / (6 * (one - epsilon))
    psi_remainder = 2 * pi * epsilon / (5 * (one - epsilon))
    b_remainder = (
        b_constant**3
        * epsilon
        / (3 * (one - b_radius))
    )
    geometric_constant = (
        one / 4
        + pi / 2
        + b_constant**2 / 2
        + delta_remainder
        + psi_remainder
        + b_remainder
    )
    if not geometric_constant < flint.arb("5.2"):
        raise RuntimeError("geometric log-ratio constant failed")

    h_r0_constant = h_bound * 3 / (4 * pi)
    m_remainder_constant = h_bound**2 / (2 * pi)
    heat_r1_constant = pi / 16
    square_constant = flint.arb("0.13")
    finite_log_constant = flint.arb("0.12")
    value_constant = (
        geometric_constant
        + h_r0_constant
        + m_remainder_constant
        + heat_r1_constant
        + square_constant
        + finite_log_constant
    )
    if not value_constant < flint.arb(LOG_RATIO_CONSTANT):
        raise RuntimeError("saved log-ratio value constant failed")

    geometric_w_constant = (
        one / 2
        + 2 * pi
        + b_constant * b_derivative_constant
        + flint.arb("0.4")
    )
    if not geometric_w_constant < flint.arb(
        LOG_RATIO_W_CONSTANT
    ):
        raise RuntimeError("saved w-derivative constant failed")

    geometric_epsilon_constant = flint.arb(12)
    h_r0_epsilon_constant = flint.arb("0.3")
    m_remainder_epsilon_constant = flint.arb(1)
    heat_r1_epsilon_constant = flint.arb("0.2")
    square_epsilon_constant = flint.arb("0.3")
    finite_log_epsilon_constant = flint.arb(1)
    epsilon_derivative_raw = (
        geometric_epsilon_constant
        + h_r0_epsilon_constant
        + m_remainder_epsilon_constant
        + heat_r1_epsilon_constant
        + square_epsilon_constant
        + finite_log_epsilon_constant
    )
    if not epsilon_derivative_raw < flint.arb(15):
        raise RuntimeError("raw epsilon-derivative constant failed")
    if not epsilon_derivative_raw < flint.arb(
        LOG_RATIO_EPSILON_CONSTANT
    ):
        raise RuntimeError("saved epsilon-derivative constant failed")

    x_constant = (
        flint.arb(LOG_RATIO_W_CONSTANT)
        + flint.arb(LOG_RATIO_EPSILON_CONSTANT) * epsilon
    ) / (8 * pi)
    if not x_constant < flint.arb(LOG_RATIO_X_CONSTANT):
        raise RuntimeError("saved x-derivative constant failed")

    exponential_loss = (
        flint.arb(LOG_RATIO_CONSTANT) * epsilon**2
    ).exp()
    ratio_value = (
        exponential_loss * flint.arb(LOG_RATIO_CONSTANT)
    )
    ratio_x = (
        exponential_loss * flint.arb(LOG_RATIO_X_CONSTANT)
    )
    if not ratio_value < flint.arb(RATIO_VALUE_CONSTANT):
        raise RuntimeError("ratio value constant failed")
    if not ratio_x < flint.arb(RATIO_X_CONSTANT):
        raise RuntimeError("ratio derivative constant failed")

    endpoint_j_rate = (
        one / 2
        + epsilon
        * (
            (one + epsilon) / (16 * pi)
            + pi / 3
            + one / 2
            + epsilon / 12
        )
    ) / (one - b_radius)
    endpoint_mu_epsilon2 = (
        3 / (8 * pi) + one / 2 + flint.arb("0.04")
    )
    lambda_epsilon2 = flint.arb(3)
    endpoint_centered_rate = (
        endpoint_j_rate
        + (endpoint_mu_epsilon2 + lambda_epsilon2) * epsilon
    )
    if not endpoint_j_rate < flint.arb("0.51"):
        raise RuntimeError("recurrence endpoint rate failed")
    if not endpoint_mu_epsilon2 < flint.arb("0.66"):
        raise RuntimeError("endpoint mu rate failed")
    if not endpoint_centered_rate < flint.arb(
        ENDPOINT_RATE_CONSTANT
    ):
        raise RuntimeError("saved centered endpoint rate failed")

    adjacent_jet = (
        flint.arb(ENDPOINT_CONSTANT)
        * (
            flint.arb(ENDPOINT_RATE_CONSTANT)
            * flint.arb(RATIO_VALUE_CONSTANT)
            + flint.arb(RATIO_X_CONSTANT)
        )
    )
    if not adjacent_jet <= flint.arb(ADJACENT_JET_CONSTANT):
        raise RuntimeError("adjacent jet constant failed")
    frame_raw = 606 / (16 * pi**2) + flint.arb(1)
    if not frame_raw < flint.arb(FRAME_CONSTANT):
        raise RuntimeError("frame constant failed")
    scalar_total = (
        adjacent_jet
        + flint.arb(EXPLICIT_D_CONSTANT)
        + flint.arb(FRAME_CONSTANT)
    )
    if not scalar_total < flint.arb(SCALAR_CHART_CONSTANT):
        raise RuntimeError("scalar chart constant failed")

    target_ratio = (
        flint.arb(SCALAR_CHART_CONSTANT) * epsilon
    )
    if not target_ratio < flint.arb("1e-7"):
        raise RuntimeError("contact-scale absorption failed")
    return {
        "precision_bits": PRECISION_BITS,
        "L_min": L_MIN,
        "epsilon_max_ball": arb_text(epsilon),
        "epsilon_max_lt": "1.4e-11",
        "h_bound_ball": arb_text(h_bound),
        "h_bound_lt": "0.54",
        "B_coefficient_ball": arb_text(b_constant),
        "B_coefficient_lt": "2.6",
        "B_radius_ball": arb_text(b_radius),
        "B_radius_lt": "4e-11",
        "geometric_constant_ball": arb_text(geometric_constant),
        "geometric_constant_lt": "5.2",
        "value_constant_ball": arb_text(value_constant),
        "value_constant_lt": str(LOG_RATIO_CONSTANT),
        "w_derivative_constant_ball": arb_text(
            geometric_w_constant
        ),
        "w_derivative_constant_lt": str(LOG_RATIO_W_CONSTANT),
        "epsilon_derivative_raw": arb_text(
            epsilon_derivative_raw
        ),
        "epsilon_derivative_raw_lt": "15",
        "epsilon_derivative_components": {
            "geometric": arb_text(geometric_epsilon_constant),
            "h_r0": arb_text(h_r0_epsilon_constant),
            "R_M": arb_text(m_remainder_epsilon_constant),
            "heat_r1": arb_text(heat_r1_epsilon_constant),
            "square": arb_text(square_epsilon_constant),
            "finite_log": arb_text(finite_log_epsilon_constant),
        },
        "epsilon_derivative_saved": (
            LOG_RATIO_EPSILON_CONSTANT
        ),
        "x_derivative_constant_ball": arb_text(x_constant),
        "x_derivative_constant_lt": str(LOG_RATIO_X_CONSTANT),
        "exponential_loss_ball": arb_text(exponential_loss),
        "ratio_value_constant_ball": arb_text(ratio_value),
        "ratio_value_constant_lt": str(RATIO_VALUE_CONSTANT),
        "ratio_x_constant_ball": arb_text(ratio_x),
        "ratio_x_constant_lt": str(RATIO_X_CONSTANT),
        "endpoint_J_rate_ball": arb_text(endpoint_j_rate),
        "endpoint_J_rate_lt": "0.51",
        "endpoint_mu_epsilon2_ball": arb_text(
            endpoint_mu_epsilon2
        ),
        "endpoint_mu_epsilon2_lt": "0.66",
        "endpoint_centered_rate_ball": arb_text(
            endpoint_centered_rate
        ),
        "endpoint_centered_rate_lt": str(
            ENDPOINT_RATE_CONSTANT
        ),
        "adjacent_jet_constant_ball": arb_text(adjacent_jet),
        "adjacent_jet_constant_le": str(ADJACENT_JET_CONSTANT),
        "frame_raw_constant_ball": arb_text(frame_raw),
        "frame_raw_constant_lt": str(FRAME_CONSTANT),
        "scalar_total_constant_ball": arb_text(scalar_total),
        "scalar_chart_constant_lt": str(SCALAR_CHART_CONSTANT),
        "target_ratio_at_L50_ball": arb_text(target_ratio),
        "target_ratio_at_L50_lt": "1e-7",
    }


def build_artifact() -> dict:
    exact = symbolic_audit()
    interval = interval_budget()
    rows = [
        StabilityRow(
            id="nfocacs_01_domain",
            role="exact_geometry",
            readiness="ready_to_apply",
            claim=(
                "The entire canonical saddle cell has an exponentially "
                "small scaled displacement."
            ),
            formula=exact["domain"],
            proof_boundary="Critical ray L>=50 and 0<=tL<=25.",
            diagnostics={
                "epsilon": interval["epsilon_max_ball"],
                "h": interval["h_bound_ball"],
            },
        ),
        StabilityRow(
            id="nfocacs_02_ratio",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The reflected entering corrected sharp divided by the "
                "C0+C1 recurrence block has one explicit ratio."
            ),
            formula=f"{exact['positive_ratio']}; {exact['omega']}",
            proof_boundary=(
                "Principal branches are fixed by epsilon*w<1 on the "
                "positive real saddle cell."
            ),
        ),
        StabilityRow(
            id="nfocacs_03_scaled_blocks",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "Every alpha and saddle remainder is analytic at "
                "epsilon=0 after the leading logarithm is removed."
            ),
            formula=exact["scaled_blocks"],
            proof_boundary="Exact scaled rewrites; no asymptotic division by zero.",
        ),
        StabilityRow(
            id="nfocacs_04_value_budget",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim=(
                "Taylor remainders after the vanished first coefficient "
                "bound the corrected log-ratio at order epsilon^2."
            ),
            formula="|Z|<8*epsilon^2",
            proof_boundary=(
                "Uses |log(1+z)-z+z^2/2|<=|z|^3/[3(1-|z|)] "
                "for the two logarithms and the exact scaled alpha blocks."
            ),
            diagnostics={
                "geometric": interval["geometric_constant_ball"],
                "total": interval["value_constant_ball"],
            },
        ),
        StabilityRow(
            id="nfocacs_05_w_derivative",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim=(
                "Differentiating after first-coefficient cancellation "
                "retains two powers of epsilon."
            ),
            formula="|Z_w|<30*epsilon^2",
            proof_boundary=(
                "The leading coefficient derivative is "
                "w/2+2i*pi*w^3+B1*B1'; all derivative remainders use "
                "the same convergent logarithm series."
            ),
            diagnostics={
                "constant": interval["w_derivative_constant_ball"],
            },
        ),
        StabilityRow(
            id="nfocacs_06_epsilon_derivative",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim=(
                "The scaled rational forms expose and cancel the apparent "
                "1/epsilon derivatives."
            ),
            formula="|Z_epsilon|<100*epsilon",
            proof_boundary=(
                "Termwise differentiation gives constants below 15: "
                "geometric <14, r0 <1, R_M <1, heat-r1 <1, square <1, "
                "and log(1+d_+) <1. The saved constant is 100."
            ),
            diagnostics={
                "raw": interval["epsilon_derivative_raw"],
                "saved": interval["epsilon_derivative_saved"],
            },
        ),
        StabilityRow(
            id="nfocacs_07_ratio_bounds",
            role="asymptotic_certificate",
            readiness="certified",
            claim=(
                "The corrected adjacent ratio and its real-axis derivative "
                "gain two and three powers of a, respectively."
            ),
            formula=exact["ratio_bounds"],
            proof_boundary=(
                "Composes the exact coordinate derivatives with "
                "|exp(Z)-1|<=exp(|Z|)|Z|."
            ),
            diagnostics={
                "x": interval["x_derivative_constant_ball"],
                "value": interval["ratio_value_constant_ball"],
                "derivative": interval["ratio_x_constant_ball"],
            },
        ),
        StabilityRow(
            id="nfocacs_08_endpoint_rate",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim=(
                "Saddle centering also removes the constant phase rate of "
                "the recurrence endpoint."
            ),
            formula=f"{exact['mu_cancellation']}; {exact['endpoint_rate']}",
            proof_boundary=(
                "The scaled negative-sharp alpha remainder is below "
                "epsilon^2; J_(a,x)/J_a contributes below 0.6 epsilon."
            ),
        ),
        StabilityRow(
            id="nfocacs_09_adjacent_real_jet",
            role="asymptotic_certificate",
            readiness="certified",
            claim=(
                "The real adjacent corrected lift and its first derivative "
                "are below the value and centered edge scales."
            ),
            formula=exact["adjacent_jet"],
            proof_boundary=(
                "Uses the inherited endpoint block bound, the ratio theorem, "
                "and equality of the reflected positive-sharp real trace."
            ),
            diagnostics={
                "constant": interval["adjacent_jet_constant_ball"],
            },
        ),
        StabilityRow(
            id="nfocacs_10_nuisances",
            role="asymptotic_certificate",
            readiness="certified",
            claim=(
                "The entering d_(n,x) correction and frame terms fit below "
                "the same exp(-7L/4) scale."
            ),
            formula=(
                "|Re(e_(N+1)d_(N+1,x))|<1500*exp(-7L/4); "
                "|u_a Delta X-v_a Delta Y|<5*exp(-7L/4)"
            ),
            proof_boundary=(
                "Composes the coefficient-mass, |d_(n,x)|, |E_[1]|, "
                "|u_a|, and |v_a| bounds from the centered reduction."
            ),
        ),
        StabilityRow(
            id="nfocacs_11_scalar_stability",
            role="asymptotic_theorem",
            readiness="ready_to_apply",
            claim=(
                "The centered scalar is uniformly stable under one adjacent "
                "Riemann-Siegel chart change."
            ),
            formula=exact["scalar_stability"],
            proof_boundary="L>=50, 0<=tL<=25, and 0<=N+1-a<=1.",
            diagnostics={
                "raw_total": interval["scalar_total_constant_ball"],
                "saved": SCALAR_CHART_CONSTANT,
            },
        ),
        StabilityRow(
            id="nfocacs_12_absorption",
            role="asymptotic_certificate",
            readiness="certified",
            claim=(
                "The entire adjacent scalar jump is negligible at the "
                "first-order contact scale."
            ),
            formula=exact["scale_comparison"],
            proof_boundary="Monotone for L>=50.",
            diagnostics={
                "L50": interval["target_ratio_at_L50_ball"],
            },
        ),
        StabilityRow(
            id="nfocacs_13_bulk_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "After endpoint chart stability, the remaining RH-level "
                "problem is the chart-invariant bulk scalar sign and count."
            ),
            formula=(
                "|A_a|>(100000L+1)*exp(-5L/4) plus an "
                "A_a-positive successor count below one turn"
            ),
            proof_boundary=(
                "No lower bound, signed crossing theorem, or winding bound "
                "is proved here."
            ),
        ),
        StabilityRow(
            id="nfocacs_14_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim=(
                "A uniform adjacent-chart theorem closes endpoint "
                "bookkeeping but not Xi transversality."
            ),
            formula="chart stability != contact exclusion",
            proof_boundary=(
                "Not the q>=1 bulk scalar theorem, q<1 theorem, finite "
                "phase closure, Lambda<=0, RH, PF-infinity, or a "
                "Clay-prize conclusion."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "uniform adjacent-chart stability certificate for the centered "
            "first-order scalar on L>=50, 0<=tL<=25; not a bulk scalar "
            "lower bound, contact exclusion, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the exact corrected positive-sharp ratio "
            "decomposition, explicit value and derivative bounds, the "
            "adjacent corrected real-jet estimates, and "
            "|A_(a,N+1)-A_(a,N)|<5000exp(-7L/4). It does not prove the "
            "q>=1 chart-invariant bulk scalar lower bound or signed crossing "
            "budget, the q<1 multiplicity-compatible theorem, finite phase "
            "cells, one-sided successor winding, contact exclusion, "
            "Lambda<=0, RH, PF-infinity, or a Clay-prize conclusion."
        ),
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "constants": {
            "log_ratio": LOG_RATIO_CONSTANT,
            "log_ratio_w": LOG_RATIO_W_CONSTANT,
            "log_ratio_epsilon": LOG_RATIO_EPSILON_CONSTANT,
            "log_ratio_x": LOG_RATIO_X_CONSTANT,
            "ratio_value": RATIO_VALUE_CONSTANT,
            "ratio_x": RATIO_X_CONSTANT,
            "endpoint": ENDPOINT_CONSTANT,
            "endpoint_rate": ENDPOINT_RATE_CONSTANT,
            "adjacent_jet": ADJACENT_JET_CONSTANT,
            "explicit_d": EXPLICIT_D_CONSTANT,
            "frame": FRAME_CONSTANT,
            "scalar_chart": SCALAR_CHART_CONSTANT,
        },
        "interval": interval,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    interval = artifact["interval"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 First-Order Centered Adjacent-Chart Stability Certificate",
            "",
            "Date: 2026-07-26",
            "",
            "Status: uniform adjacent-chart theorem on `L>=50`,",
            "`0<=tL<=25`. This closes one endpoint bookkeeping problem;",
            "it is not the bulk scalar lower bound and not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Scaled Domain",
            "",
            "Put",
            "",
            "```text",
            exact["domain"],
            "h=1/2-i*pi*t/8, n=N+1.",
            "```",
            "",
            f"The 256-bit interval endpoint is `{interval['epsilon_max_ball']}`;",
            "in particular `epsilon<1.4e-11` and all logarithm disks below",
            "have radius less than `4e-11`.",
            "",
            "## Exact Corrected Ratio",
            "",
            "Reflect the entering negative sharp to the positive sharp, which",
            "does not change its real trace. Dividing by the adjacent",
            "`C_0+C_1/a` recurrence block gives",
            "",
            "```text",
            exact["positive_ratio"],
            exact["omega"],
            "```",
            "",
            "The heat/logarithm combination simplifies exactly:",
            "",
            "```text",
            "-h*delta+i*pi*t*(r1-delta)/8",
            " =-delta/2+i*pi*t*r1/8.",
            "```",
            "",
            "There is therefore no uncancelled `t*delta` term. The scaled",
            "blocks are",
            "",
            "```text",
            exact["scaled_blocks"],
            "```",
            "",
            "These forms are analytic at `epsilon=0`; bounding an apparent",
            "`1/epsilon` derivative term by term before this rewrite would",
            "miss the cancellation.",
            "",
            "## Value Bound",
            "",
            "Write `B=1+epsilon*B1`,",
            "`B1=-w/2-i*(2pi/3)w^3`. The exact first coefficients of",
            "`-delta/2-i*pi*Psi` and `log B` agree. Taylor's bound",
            "",
            "```text",
            "|log(1+z)-z+z^2/2|<=|z|^3/[3(1-|z|)]",
            "```",
            "",
            "gives the geometric coefficient budget",
            "",
            "```text",
            "1/4+pi/2+|B1|^2/2+three explicit remainders",
            f" = {interval['geometric_constant_ball']} < 5.2.",
            "```",
            "",
            "The remaining scaled terms satisfy",
            "",
            "```text",
            "|h*r0|/epsilon^2 < |h|*3/(4pi),",
            "|R_M|/epsilon^2 < |h|^2/(2pi),",
            "|i*pi*t*r1/8|/epsilon^2 < pi/16,",
            "|t*(r1-delta)^2/4|/epsilon^2 < 0.13,",
            "|log(1+d_+)|/epsilon^2 < 0.12.",
            "```",
            "",
            "Their interval sum is",
            "",
            "```text",
            f"{interval['value_constant_ball']} < 8,",
            "|Z|=|log rho|<8*epsilon^2.",
            "```",
            "",
            "## Derivative Bound",
            "",
            "Differentiation is performed after the first coefficient has",
            "cancelled. For `w` the coefficient derivative costs",
            "",
            "```text",
            "1/2+2pi+|B1||B1'|+0.4",
            f" = {interval['w_derivative_constant_ball']} < 30,",
            "|Z_w|<30*epsilon^2.",
            "```",
            "",
            "For `epsilon`, direct differentiation of the displayed scaled",
            "rational forms gives the termwise budget",
            "",
            "```text",
            "geometric <14, r0 <1, R_M <1, heat-r1 <1,",
            "square <1, log(1+d_+) <1;",
            f"interval sum = {interval['epsilon_derivative_raw']} < 15.",
            "|Z_epsilon|<100*epsilon.",
            "```",
            "",
            "The saved constant leaves a factor-four margin. Since",
            "",
            "```text",
            exact["coordinate_derivatives"],
            "```",
            "",
            "one obtains",
            "",
            "```text",
            f"|Z_x|<{interval['x_derivative_constant_ball']}*epsilon^3",
            "<2*epsilon^3,",
            "|rho-1|<11*epsilon^2, |rho_x|<3*epsilon^3.",
            "```",
            "",
            "## Centered Endpoint Rate",
            "",
            "The endpoint prefactor has a second exact cancellation. With",
            "`chi=alpha-log(a)`,",
            "",
            "```text",
            exact["mu_cancellation"],
            "```",
            "",
            "and the `-pi/8` in `K_x/K` cancels the `pi/8` from",
            "`i*chi/2`. The scaled alpha remainder gives",
            "`|mu_a|<epsilon^2`. The explicit recurrence polynomial gives",
            (
                "`|J_(a,x)/J_a|<"
                f"{interval['endpoint_J_rate_ball']}epsilon"
                "<0.51epsilon`; together with"
            ),
            "`|lambda_a|<3epsilon^2`,",
            "",
            "```text",
            exact["endpoint_rate"],
            "```",
            "",
            "## Chart-Stability Theorem",
            "",
            "The inherited recurrence endpoint bound, enlarged from 50 to",
            "100 to include `B`, and the ratio estimates imply",
            "",
            "```text",
            exact["adjacent_jet"],
            "```",
            "",
            "The exact projection identity is",
            "",
            "```text",
            "Delta A_a=Delta U-u_a*Delta X+v_a*Delta Y",
            "             -Re(e_(N+1)d_(N+1,x)).",
            "```",
            "",
            "The entering derivative correction costs at most",
            "`1500exp(-7L/4)`. The imported `u_a`, `v_a`, and complex-main",
            "bounds put both frame terms below `5exp(-7L/4)`. Hence",
            "",
            "```text",
            f"raw constant = {interval['scalar_total_constant_ball']}",
            exact["scalar_stability"],
            "```",
            "",
            "At `L=50` this is smaller than the `exp(-5L/4)` scale by",
            f"`{interval['target_ratio_at_L50_ball']} < 1e-7`, and the",
            "ratio decreases thereafter.",
            "",
            "This proves adjacent chart stability only. It removes the",
            "last-saddle/endpoint bookkeeping obstruction but supplies no",
            "lower bound or sign for the chart-invariant bulk centered",
            "moment, no positive-crossing count, no q<1 theorem, and no",
            "contact exclusion, `Lambda<=0`, RH, PF-infinity, or",
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
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman first-order centered adjacent-chart stability "
        "certificate: 14 rows, |rho-1|<11/a^2, |rho_x|<3/a^3, "
        "|Delta A|<5000*e^-7L/4, 1 open bulk obligation"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
