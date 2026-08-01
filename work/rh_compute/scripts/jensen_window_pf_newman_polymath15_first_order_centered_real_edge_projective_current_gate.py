#!/usr/bin/env python3
"""Build the q=1 real-edge projective-current asymptotic gate."""

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
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_real_edge_projective_current_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "endpoint_relative_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_endpoint_relative_phase_current_"
        "recurrence_gate.json"
    ),
    "complex_endpoint_source": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complex_endpoint_source_"
        "normalization_gate.json"
    ),
    "adjacent_chart_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_chart_stability_certificate.json"
    ),
    "interior_projective_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_interior_projective_current_gate.json"
    ),
    "c0_source_extraction": (
        REPO_ROOT
        / "work/rh_compute/scripts/"
        "jensen_window_pf_newman_polymath15_"
        "critical_RS_C1_endpoint_peeling_contract.py"
    ),
}

PRECISION_BITS = 192
SINC_TERMS = 24
CURVATURE_MARGIN = Fraction(3, 50)
CURRENT_MARGIN = Fraction(3, 8000)
REGIONS = (
    (Fraction(-1), Fraction(-5, 8), 512, None),
    (Fraction(-5, 8), Fraction(-3, 8), 512, -1),
    (Fraction(-3, 8), Fraction(3, 8), 1024, None),
    (Fraction(3, 8), Fraction(5, 8), 512, 1),
    (Fraction(5, 8), Fraction(1), 512, None),
)


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
            "location": "equation (53), C_0(p), and Proposition 6.2",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def aq(value: Fraction | int) -> flint.arb:
    value = Fraction(value)
    return flint.arb(flint.fmpq(value.numerator, value.denominator))


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def mul_jet(
    left: tuple[flint.arb, flint.arb, flint.arb],
    right: tuple[flint.arb, flint.arb, flint.arb],
) -> tuple[flint.arb, flint.arb, flint.arb]:
    a0, a1, a2 = left
    b0, b1, b2 = right
    return (
        a0 * b0,
        a1 * b0 + a0 * b1,
        a2 * b0 + 2 * a1 * b1 + a0 * b2,
    )


def inv_jet(
    value: tuple[flint.arb, flint.arb, flint.arb],
) -> tuple[flint.arb, flint.arb, flint.arb]:
    f0, f1, f2 = value
    return (
        1 / f0,
        -f1 / f0**2,
        2 * f1**2 / f0**3 - f2 / f0**2,
    )


def sinc_derivative_balls(
    z: flint.arb,
) -> tuple[flint.arb, flint.arb, flint.arb]:
    """Enclose sinc(z) and its first two derivatives by entire series."""
    upper = flint.arb(abs(z).upper())
    value = flint.arb(0)
    first = flint.arb(0)
    second = flint.arb(0)
    for k in range(SINC_TERMS):
        sign = -1 if k % 2 else 1
        value += sign * z ** (2 * k) / factorial(2 * k + 1)
    for k in range(1, SINC_TERMS):
        sign = -1 if k % 2 else 1
        first += (
            sign * (2 * k) * z ** (2 * k - 1) / factorial(2 * k + 1)
        )
        second += (
            sign
            * (2 * k)
            * (2 * k - 1)
            * z ** (2 * k - 2)
            / factorial(2 * k + 1)
        )

    k = SINC_TERMS
    ratio_0 = upper**2 / ((2 * k + 2) * (2 * k + 3))
    ratio_1 = upper**2 / (2 * k * (2 * k + 3))
    ratio_2 = (
        (2 * k + 1)
        * upper**2
        / (2 * k * (2 * k - 1) * (2 * k + 3))
    )
    if not ratio_0 < 1 or not ratio_1 < 1 or not ratio_2 < 1:
        raise RuntimeError("sinc tail ratio is not contractive")
    error_0 = (
        upper ** (2 * k)
        / factorial(2 * k + 1)
        / (1 - ratio_0)
    )
    error_1 = (
        (2 * k)
        * upper ** (2 * k - 1)
        / factorial(2 * k + 1)
        / (1 - ratio_1)
    )
    error_2 = (
        (2 * k)
        * (2 * k - 1)
        * upper ** (2 * k - 2)
        / factorial(2 * k + 1)
        / (1 - ratio_2)
    )
    return (
        value + flint.arb(0, error_0),
        first + flint.arb(0, error_1),
        second + flint.arb(0, error_2),
    )


def sinc_composed_jet(
    z: tuple[flint.arb, flint.arb, flint.arb],
) -> tuple[flint.arb, flint.arb, flint.arb]:
    z0, z1, z2 = z
    s0, s1, s2 = sinc_derivative_balls(z0)
    return (s0, s1 * z1, s2 * z1**2 + s1 * z2)


def local_a_jet(
    p: flint.arb,
    center_sign: int,
) -> tuple[flint.arb, flint.arb, flint.arb]:
    """Evaluate the removable A jet in a sinc chart at p=+-1/2."""
    pi = flint.arb.pi()
    y = p - aq(Fraction(center_sign, 2))
    if center_sign == 1:
        delta = (pi / 2) * y * (y - 3)
        linear = ((y - 3) / 4, flint.arb(1) / 4, flint.arb(0))
        delta_jet = (delta, (pi / 2) * (2 * y - 3), pi)
    elif center_sign == -1:
        delta = (pi / 2) * y * (y - 5)
        linear = ((5 - y) / 4, -flint.arb(1) / 4, flint.arb(0))
        delta_jet = (delta, (pi / 2) * (2 * y - 5), pi)
    else:
        raise ValueError("center_sign must be -1 or 1")
    denominator_jet = (pi * y, pi, flint.arb(0))
    quotient = mul_jet(
        sinc_composed_jet(delta_jet),
        inv_jet(sinc_composed_jet(denominator_jet)),
    )
    return mul_jet(linear, quotient)


def direct_a_jet(
    p: flint.arb,
) -> tuple[flint.arb, flint.arb, flint.arb]:
    variable = flint.arb_series([p, flint.arb(1)], prec=3)
    phase = variable * variable / 2 - 2 * variable + flint.arb(3) / 8
    value = -phase.cos_pi() / (2 * variable.cos_pi())
    coefficients = value.coeffs()
    return coefficients[0], coefficients[1], 2 * coefficients[2]


def interval_ball(left: Fraction, right: Fraction) -> flint.arb:
    lower = aq(left)
    upper = aq(right)
    return (lower + upper) / 2 + flint.arb(0, (upper - lower) / 2)


def point_diagnostics() -> dict:
    pi = flint.arb.pi()
    rows = {}
    for label, point in (
        ("midpoint", Fraction(0)),
        ("phase_guard", Fraction(29, 30)),
        ("cutoff", Fraction(1)),
    ):
        a0, a1, a2 = direct_a_jet(aq(point))
        curvature = a1**2 - a0 * a2
        current = -curvature / (16 * pi**2)
        rows[label] = {
            "p": str(point),
            "A_ball": arb_text(a0),
            "curvature_ball": arb_text(curvature),
            "current_symbol_ball": arb_text(current),
        }

    p = aq(Fraction(29, 30))
    variable = flint.arb_series([p, flint.arb(1)], prec=3)
    psi = variable * variable / 2 + flint.arb(3) / 8
    real_c0 = psi.cos_pi() / (2 * variable.cos_pi())
    f0, f1, f2_half = real_c0.coeffs()
    phi = p * p / 2 + p + flint.arb(3) / 8
    theta = (1 - p) / 2
    wrong_a = f0 - phi.cos_pi()
    wrong_b = -f1 / (4 * pi) + theta * phi.sin_pi() / 2
    wrong_m = (
        2 * f2_half / (16 * pi**2)
        + phi.sin_pi() / (16 * pi)
        + theta**2 * phi.cos_pi() / 4
    )
    wrong_current = wrong_a * wrong_m - wrong_b**2
    rows["phase_guard"]["wrong_terminal_phase_current_ball"] = arb_text(
        wrong_current
    )
    return rows


def interval_certificate() -> dict:
    flint.ctx.prec = PRECISION_BITS
    margin = aq(CURVATURE_MARGIN)
    region_rows = []
    global_lower: flint.arb | None = None
    global_location = ""
    total = 0
    for left, right, boxes, center_sign in REGIONS:
        width = (right - left) / boxes
        region_lower: flint.arb | None = None
        region_location = ""
        for index in range(boxes):
            box_left = left + index * width
            box_right = box_left + width
            p = interval_ball(box_left, box_right)
            jet = (
                local_a_jet(p, center_sign)
                if center_sign is not None
                else direct_a_jet(p)
            )
            curvature = jet[1] ** 2 - jet[0] * jet[2]
            if not curvature > margin:
                raise RuntimeError(
                    f"uncertified curvature box {box_left}..{box_right}: "
                    f"{curvature}"
                )
            lower_ball = flint.arb(curvature.lower())
            if region_lower is None or lower_ball < region_lower:
                region_lower = lower_ball
                region_location = f"{box_left}..{box_right}"
            if global_lower is None or lower_ball < global_lower:
                global_lower = lower_ball
                global_location = f"{box_left}..{box_right}"
            total += 1
        assert region_lower is not None
        region_rows.append(
            {
                "interval": f"[{left},{right}]",
                "chart": (
                    "direct quotient"
                    if center_sign is None
                    else f"sinc chart at {center_sign}/2"
                ),
                "boxes": boxes,
                "minimum_lower_ball": arb_text(region_lower),
                "minimum_box": region_location,
            }
        )

    pi = flint.arb.pi()
    if not pi**2 < 10:
        raise RuntimeError("pi^2<10 audit failed")
    assert global_lower is not None
    return {
        "precision_bits": PRECISION_BITS,
        "sinc_terms": SINC_TERMS,
        "domain": "-1<=p<=1",
        "regions": region_rows,
        "boxes": total,
        "removable_centers": ["-1/2", "1/2"],
        "curvature": "D_edge^(0)=(A')^2-A*A''",
        "curvature_margin": ">3/50",
        "minimum_lower_ball": arb_text(global_lower),
        "minimum_box": global_location,
        "pi_square_bound": "pi^2<10",
        "current_symbol_margin": "K_edge<-3/8000",
        "points": point_diagnostics(),
    }


def verify_symbolics() -> None:
    c, cx, gr, grx, alpha, alpha_x = sp.symbols(
        "C C_x G_R G_xR alpha alpha_x", real=True
    )
    d = gr - alpha * c
    dx = grx - alpha_x * c - alpha * cx
    current = sp.expand(c * dx - d * cx)
    target = sp.expand(c * grx - gr * cx - alpha_x * c**2)
    if sp.simplify(current - target) != 0:
        raise RuntimeError("division-free edge-current cancellation failed")

    p = sp.symbols("p", real=True)
    pi = sp.pi
    theta = (1 - p) / 2
    psi = pi * (p**2 / 2 + sp.Rational(3, 8))
    delta = pi * (p**2 / 2 - p + sp.Rational(3, 8))
    f = sp.cos(psi) / (2 * sp.cos(pi * p))
    a = f - sp.cos(delta)
    b = -sp.diff(f, p) / (4 * pi) + theta * sp.sin(delta) / 2
    m = (
        sp.diff(f, p, 2) / (16 * pi**2)
        + sp.sin(delta) / (16 * pi)
        + theta**2 * sp.cos(delta) / 4
    )
    if sp.simplify(b + sp.diff(a, p) / (4 * pi)) != 0:
        raise RuntimeError("first edge-jet collapse failed")
    if sp.simplify(m - sp.diff(a, p, 2) / (16 * pi**2)) != 0:
        raise RuntimeError("second edge-jet collapse failed")
    current_symbol = sp.simplify(a * m - b**2)
    curvature_symbol = sp.simplify(
        (a * sp.diff(a, p, 2) - sp.diff(a, p) ** 2)
        / (16 * pi**2)
    )
    if sp.simplify(current_symbol - curvature_symbol) != 0:
        raise RuntimeError("edge log-curvature identity failed")

    h = psi - 2 * pi * p
    shifted = -sp.cos(h) / (2 * sp.cos(pi * p))
    if sp.trigsimp(a - shifted, method="fu") != 0:
        raise RuntimeError("shifted real C_0 identity failed")

    scale, radial, aa, bb, mm = sp.symbols(
        "S r_S A B M", real=True, nonzero=True
    )
    edge_c = scale * aa
    edge_d = scale * bb / sp.Symbol("a", positive=True)
    saddle = next(symbol for symbol in edge_d.free_symbols if symbol.name == "a")
    edge_cx = radial * edge_c + scale * bb / saddle
    edge_dx = radial * edge_d + scale * mm / saddle**2
    normalized = sp.simplify(
        saddle**2 * (edge_c * edge_dx - edge_d * edge_cx) / scale**2
    )
    if sp.simplify(normalized - (aa * mm - bb**2)) != 0:
        raise RuntimeError("radial edge-current limit failed")


def exact_payload() -> dict:
    return {
        "edge_current": {
            "coordinates": (
                "V=e+z_N, G=g+s_*'*u_N*z_N, C_edge=Re(V), "
                "D_edge=Re(G)-alpha*C_edge, alpha=Re(s_*')*u_N"
            ),
            "division_free": (
                "J_edge=C_edge*D_(edge,x)-D_edge*C_(edge,x)="
                "C_edge*Re(G_x)-Re(G)*C_(edge,x)-alpha_x*C_edge^2"
            ),
            "retained_second_jet": (
                "G_x=g_x+(s_*''u_N+s_*'u_(N,x)+"
                "s_*'u_N*z_(N,x)/z_N)z_N; "
                "g_x retains H_(a,xx), mu_(a,x), the moving T_0+i "
                "factor, and the real radial kappa rate"
            ),
        },
        "q1_source": {
            "path": (
                "theta=(1-p)/2, a=N+theta, q=2tL^2=1, "
                "S_a=kappa*T_0, r_S=S_(a,x)/S_a, fixed p in [-1,1], "
                "N->infinity"
            ),
            "endpoint": (
                "F=C_0(p), H_a->F, a*H_(a,x)->-F'/(4pi), "
                "a^2*H_(a,xx)->F''/(16pi^2), a*mu_a->0, "
                "a^2*mu_(a,x)->0"
            ),
            "terminal_phase": (
                "Q(p)=exp[-pi*i*(p^2/2-p+3/8)], "
                "z_N/S_a->-Q(p), a*u_N->theta, "
                "a*(z_(N,x)/z_N-r_S)->-i*theta/2"
            ),
            "phase_guard": (
                "Q=conjugate(R)*exp[-2pi*i*(1-p)] for "
                "R=exp[pi*i*(p^2/2+p+3/8)]. The extra factor is one at "
                "p=0 and p=1 but not in the interior; endpoint-only tests "
                "cannot determine the terminal phase."
            ),
        },
        "asymptotic_jets": {
            "definitions": (
                "A=Re(F-Q); B=Re[-F'/(4pi)+i*theta*Q/2]; "
                "M=Re[F''/(16pi^2)+Q{i/(16pi)+theta^2/4}]"
            ),
            "limits": (
                "C_edge/S_a->A, a*D_edge/S_a->B, "
                "a(C_(edge,x)-r_S*C_edge)/S_a->B, and "
                "a^2(D_(edge,x)-r_S*D_edge)/S_a->M"
            ),
            "derivative_collapse": (
                "B=-A'/(4pi), M=A''/(16pi^2)"
            ),
            "current_limit": (
                "a^2*J_edge/S_a^2->K_edge(p)="
                "[A*A''-(A')^2]/(16pi^2)"
            ),
            "division_guard": (
                "The polynomial limit is primary. At A=0 it equals "
                "-(A')^2/(16pi^2), with no division by C_edge or A."
            ),
        },
        "curvature_certificate": {
            "real_source": (
                "A(p)=-cos(pi*(p^2/2-2p+3/8))/(2cos(pi*p))="
                "-Re C_0(p-2), with removable values "
                "A(-1/2)=5/4 and A(1/2)=-3/4"
            ),
            "removable_charts": (
                "At p=1/2+y, A=(y-3)sinc(pi*y*(y-3)/2)/"
                "[4sinc(pi*y)]; at p=-1/2+y, "
                "A=(5-y)sinc(pi*y*(y-5)/2)/[4sinc(pi*y)]."
            ),
            "strict_margin": (
                "(A')^2-A*A''>3/50 on -1<=p<=1; since pi^2<10, "
                "K_edge(p)<-3/8000."
            ),
        },
        "route_decision": {
            "proved_symbol": (
                "The complete q=1 cofinal real-edge leading symbol is "
                "strictly clockwise on the full saddle cell, including "
                "both removable C_0 points and every A=0 fibre."
            ),
            "next_target": (
                "Derive an explicit uniform finite-a remainder smaller "
                "than 3/8000 for a^2*J_edge/(kappa*T_0)^2 on q=1, "
                "then extend the signed edge estimate to q>=1 and splice "
                "it through adjacent cutoffs before returning to the "
                "cumulative contact scalar."
            ),
        },
    }


def build_rows(exact: dict, interval: dict) -> list[GateRow]:
    return [
        GateRow(
            "repc_01_edge_coordinates",
            "exact_definition",
            "ready_to_apply",
            "The endpoint and terminal carrier are composed before projection.",
            exact["edge_current"]["coordinates"],
            "One fixed-N chart; no sign is asserted here.",
        ),
        GateRow(
            "repc_02_division_free_current",
            "exact_identity",
            "ready_to_apply",
            "Terminal centering cancels the alpha*C*C_x terms exactly.",
            exact["edge_current"]["division_free"],
            "No division by C_edge, c_0, or H_a is used.",
        ),
        GateRow(
            "repc_03_second_jet",
            "source_identity",
            "ready_to_apply",
            "The real-edge current retains the complete endpoint and terminal second jet.",
            exact["edge_current"]["retained_second_jet"],
            "Dropping any displayed derivative changes the leading symbol.",
        ),
        GateRow(
            "repc_04_q1_terminal_phase",
            "source_asymptotic",
            "available_exact",
            "The terminal saddle has an interior phase not visible at p=0 or p=1.",
            exact["q1_source"]["terminal_phase"],
            "Fixed p on the stated q=1 cofinal family.",
        ),
        GateRow(
            "repc_05_phase_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Endpoint-only phase interpolation gives a false interior sign candidate.",
            exact["q1_source"]["phase_guard"],
            "The corrected Q phase is used in every promoted formula.",
            diagnostics=interval["points"]["phase_guard"],
        ),
        GateRow(
            "repc_06_asymptotic_jets",
            "source_asymptotic",
            "available_exact",
            "Four radial-subtracted edge jets close on three real functions.",
            exact["asymptotic_jets"]["limits"],
            "The finite-a uniform error is not yet quantified.",
        ),
        GateRow(
            "repc_07_derivative_collapse",
            "exact_identity",
            "ready_to_apply",
            "The source jets are consecutive derivatives of one real endpoint shape.",
            exact["asymptotic_jets"]["derivative_collapse"],
            "Uses the corrected terminal phase Q.",
        ),
        GateRow(
            "repc_08_current_limit",
            "asymptotic_identity",
            "available_exact",
            "The radial rate cancels and leaves one strict log-curvature symbol.",
            exact["asymptotic_jets"]["current_limit"],
            "Pointwise cofinal limit; no finite-height threshold is claimed.",
        ),
        GateRow(
            "repc_09_zero_fibre",
            "division_guard",
            "ready_to_apply",
            "The leading current remains meaningful at every zero real projection.",
            exact["asymptotic_jets"]["division_guard"],
            "No projective tangent is introduced.",
        ),
        GateRow(
            "repc_10_real_source",
            "exact_source_reduction",
            "ready_to_apply",
            "The edge shape is a shifted real C_0 trace with two removable points.",
            exact["curvature_certificate"]["real_source"],
            "The quotient is interpreted by analytic continuation.",
        ),
        GateRow(
            "repc_11_removable_charts",
            "analytic_desingularization",
            "ready_to_apply",
            "Two sinc charts remove every interval division by cos(pi*p).",
            exact["curvature_certificate"]["removable_charts"],
            "The sinc tails are enclosed as entire series.",
        ),
        GateRow(
            "repc_12_interval_cover",
            "interval_certificate",
            "certified",
            "A rational Arb cover proves strict curvature on the whole cell.",
            exact["curvature_certificate"]["strict_margin"],
            "Complete interval cover, not a floating grid.",
            diagnostics={
                "boxes": interval["boxes"],
                "precision_bits": interval["precision_bits"],
                "minimum_lower_ball": interval["minimum_lower_ball"],
                "minimum_box": interval["minimum_box"],
            },
        ),
        GateRow(
            "repc_13_uniform_symbol_sign",
            "proved_asymptotic_symbol",
            "ready_to_apply",
            "The complete q=1 cofinal edge symbol has one strict sign.",
            exact["route_decision"]["proved_symbol"],
            "A uniform finite-a remainder and q>1 extension remain open.",
        ),
        GateRow(
            "repc_14_handoff",
            "route_decision",
            "open_quantitative_handoff",
            "The leading-symbol wall is closed and the finite remainder is isolated.",
            exact["route_decision"]["next_target"],
            "No Abel gap, winding cap, contact exclusion, or RH conclusion follows yet.",
        ),
    ]


def build_payload() -> dict:
    verify_symbolics()
    exact = exact_payload()
    interval = interval_certificate()
    rows = build_rows(exact, interval)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact division-free real-edge current and rigorous q=1 "
            "cofinal negative leading-symbol certificate"
        ),
        "proof_boundary": (
            "This artifact proves the exact division-free edge-current "
            "identity, the corrected q=1 terminal phase, the fixed-p "
            "cofinal edge-current limit, and a complete Arb certificate "
            "that its leading symbol is below -3/8000 on -1<=p<=1. It "
            "does not yet prove a uniform finite-a edge sign, a q>1 edge "
            "sign, adjacent-cutoff splicing, a signed cumulative-minor "
            "estimate, an Abel-scalar gap, a successor winding cap, "
            "contact exclusion, Q209, the cofinal descendant theorem, "
            "Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
        "source_audit": source_audit(),
        "exact": exact,
        "interval": interval,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "division_free_edge_currents": 1,
            "asymptotic_edge_jet_limits": 4,
            "arb_intervals": interval["boxes"],
            "removable_charts": 2,
            "strict_curvature_margins": 1,
            "uniform_q1_leading_symbol_signs": 1,
            "uniform_finite_height_edge_signs": 0,
            "abel_gaps": 0,
            "winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    interval = payload["interval"]
    return f"""# Jensen-Window PF Newman Polymath-15 Real-Edge Projective-Current Gate

Date: 2026-07-31

Status: exact division-free edge current and rigorous negative `q=1`
cofinal leading symbol. This is not a proof of a finite-height edge theorem,
an Abel gap, `Lambda <= 0`, or RH.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Exact Edge Current

```text
{exact["edge_current"]["coordinates"]}

{exact["edge_current"]["division_free"]}
```

The `alpha*C*C_x` terms cancel. The exact derivative retains

```text
{exact["edge_current"]["retained_second_jet"]}
```

No ratio is taken at `C_edge=0`, `c_0=0`, or `H_a=0`.

## Correct Terminal Phase

On

```text
{exact["q1_source"]["path"]}
```

the endpoint and terminal limits are

```text
{exact["q1_source"]["endpoint"]}

{exact["q1_source"]["terminal_phase"]}
```

The phase guard is essential:

```text
{exact["q1_source"]["phase_guard"]}
```

At `p=29/30`, the wrong endpoint-interpolated phase gives the Arb current
ball

```text
{interval["points"]["phase_guard"]["wrong_terminal_phase_current_ball"]}
```

whereas the corrected source gives

```text
{interval["points"]["phase_guard"]["current_symbol_ball"]}
```

The full finite source was used to reject the wrong branch before promotion.

## Three-Jet Collapse

Define

```text
{exact["asymptotic_jets"]["definitions"]}
```

Then

```text
{exact["asymptotic_jets"]["limits"]}

{exact["asymptotic_jets"]["derivative_collapse"]}

{exact["asymptotic_jets"]["current_limit"]}
```

This is division-free. In particular,

```text
{exact["asymptotic_jets"]["division_guard"]}
```

## Whole-Cell Certificate

The real source is

```text
{exact["curvature_certificate"]["real_source"]}
```

At the two apparent poles the checker uses

```text
{exact["curvature_certificate"]["removable_charts"]}
```

The complete Arb cover has `{interval["boxes"]}` rational intervals at
`{interval["precision_bits"]}` bits. Its weakest saved lower endpoint is

```text
{interval["minimum_lower_ball"]}
```

on `{interval["minimum_box"]}`. Therefore

```text
{exact["curvature_certificate"]["strict_margin"]}
```

This includes every zero of `A`; no tangent chart is needed.

## Handoff

```text
{exact["route_decision"]["proved_symbol"]}

{exact["route_decision"]["next_target"]}
```

The missing finite-`a` error bound is now quantitative and has a fixed
`3/8000` reserve. The Abel gap, winding cap, contact exclusion, and RH
remain open.
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    counts = payload["counts"]
    print(
        "built real-edge projective-current gate: "
        f"{counts['rows']} rows, "
        f"{counts['division_free_edge_currents']} division-free edge current, "
        f"{counts['asymptotic_edge_jet_limits']} asymptotic edge jets, "
        f"{counts['arb_intervals']} Arb intervals, "
        f"{counts['removable_charts']} removable charts, "
        f"{counts['strict_curvature_margins']} strict curvature margin, "
        f"{counts['uniform_q1_leading_symbol_signs']} uniform q=1 leading sign, "
        f"{counts['uniform_finite_height_edge_signs']} finite-height edge signs, "
        f"{counts['abel_gaps']} Abel gaps, "
        f"{counts['winding_bounds']} winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
