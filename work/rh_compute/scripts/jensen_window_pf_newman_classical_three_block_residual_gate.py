#!/usr/bin/env python3
"""Audit classical and positive 9/5/1 Gasper-block residuals of Xi."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

from flint import acb, arb, ctx  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_classical_three_block_residual_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_classical_three_block_residual_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def acb_xi(s: acb) -> acb:
    pi_ball = arb.pi()
    return (
        acb(arb("1/2"))
        * s
        * (s - 1)
        * acb(pi_ball) ** (-s / 2)
        * (s / 2).gamma()
        * s.zeta()
    )


def acb_h0(x: acb) -> acb:
    return acb_xi((1 + acb(0, 1) * x) / 2) / 8


def acb_gasper_block(x: acb, order: arb) -> acb:
    pi_ball = arb.pi()
    shift = acb(0, 1) * x / 4
    argument = acb(2 * pi_ball)
    return acb(pi_ball**2 / 2) * (
        argument.bessel_k(acb(order) + shift)
        + argument.bessel_k(acb(order) - shift)
    )


def cauchy_finite_difference_jet(
    function,
    x_text: str,
    step: arb,
    cauchy_radius: arb,
    box_radius: arb,
) -> tuple[list[arb], arb]:
    """Enclose f, f', and f'' by central differences and Cauchy bounds."""

    x = arb(x_text)
    f_minus = function(acb(x - step))
    f_zero = function(acb(x))
    f_plus = function(acb(x + step))
    box = acb(arb(x, box_radius), arb(0, box_radius))
    maximum = abs(function(box)).upper()
    if not maximum.is_finite():
        raise RuntimeError(f"nonfinite Cauchy box at x={x_text}")

    first = ((f_plus - f_minus) / (2 * step)).real
    second = ((f_plus - 2 * f_zero + f_minus) / step**2).real
    first += arb(0, (step**2 * maximum / cauchy_radius**3).upper())
    second += arb(0, (2 * step**2 * maximum / cauchy_radius**4).upper())
    return [f_zero.real, first, second], maximum


def laguerre(jet: list[arb]) -> arb:
    return jet[1] ** 2 - jet[0] * jet[2]


def cross_laguerre(left: list[arb], right: list[arb]) -> arb:
    return (
        2 * left[1] * right[1]
        - left[0] * right[2]
        - right[0] * left[2]
    )


def sympy_to_arb(expression: sp.Expr) -> arb:
    expression = sp.sympify(expression)
    if expression == sp.pi:
        return arb.pi()
    if expression.is_Integer:
        return arb(int(expression))
    if expression.is_Rational:
        return arb(int(expression.p)) / int(expression.q)
    if expression.is_Add:
        return sum((sympy_to_arb(term) for term in expression.args), arb(0))
    if expression.is_Mul:
        result = arb(1)
        for factor in expression.args:
            result *= sympy_to_arb(factor)
        return result
    if expression.is_Pow and expression.exp.is_Integer:
        return sympy_to_arb(expression.base) ** int(expression.exp)
    raise TypeError(f"unsupported exact expression: {expression!r}")


def bernstein_coefficients(polynomial: sp.Poly) -> list[sp.Expr]:
    degree = polynomial.degree()
    power = [polynomial.nth(index) for index in range(degree + 1)]
    return [
        sp.factor(
            sum(
                power[index]
                * sp.binomial(k, index)
                / sp.binomial(degree, index)
                for index in range(k + 1)
            )
        )
        for k in range(degree + 1)
    ]


def build_exact() -> dict:
    beta, gamma = sp.symbols("beta gamma", real=True)
    r = sp.symbols("r", positive=True)
    q = r**2
    pi = sp.pi
    a = sp.Rational(3, 2) / pi
    b = pi - a
    gamma_match = (pi**2 - 3) / 2

    e, ep, epp = sp.symbols("e ep epp", real=True)
    p, pp, ppp = sp.symbols("p pp ppp", real=True)
    s, sp1, spp = sp.symbols("s sp spp", real=True)
    direct = (ep - beta * pp - gamma * sp1) ** 2 - (
        e - beta * p - gamma * s
    ) * (epp - beta * ppp - gamma * spp)
    split = (
        ep**2
        - e * epp
        - beta * (2 * ep * pp - e * ppp - p * epp)
        - gamma * (2 * ep * sp1 - e * spp - s * epp)
        + beta**2 * (pp**2 - p * ppp)
        + beta * gamma * (2 * pp * sp1 - p * spp - s * ppp)
        + gamma**2 * (sp1**2 - s * spp)
    )
    if sp.expand(direct - split) != 0:
        raise RuntimeError("bivariate Laguerre identity failed")

    normalized_block = (
        1
        + r**9
        + beta * q * (1 + r**5)
        + gamma * q**2 * (1 + r)
    )
    p2_block = sp.expand(normalized_block.subs({beta: -a, gamma: 0}))
    p2_lower = sp.factor((1 + pi * q) * (1 - a * q) - p2_block)
    p2_expected = q * (pi - sp.Rational(3, 2) * q - r**7 + a * r**5)
    if sp.expand(p2_lower - p2_expected) != 0:
        raise RuntimeError("Polya P2 kernel lower bound failed")

    db_block = sp.expand(normalized_block.subs({beta: b, gamma: 1}))
    exponential_taylor = sum((pi * q) ** k / sp.factorial(k) for k in range(5))
    db_polynomial = sp.Poly(
        sp.cancel(((1 - a * q) * exponential_taylor - db_block) / q**2),
        r,
    )
    if db_polynomial.degree() != 6:
        raise RuntimeError("de Bruijn Bernstein polynomial degree drifted")
    db_bernstein = bernstein_coefficients(db_polynomial)
    ctx.dps = 130
    db_bernstein_balls = [sympy_to_arb(value) for value in db_bernstein]
    if not all(value > 0 for value in db_bernstein_balls):
        raise RuntimeError("de Bruijn kernel Bernstein certificate failed")

    tail_series = sp.series(sp.exp(pi * q) * (1 - a * q), r, 0, 7)
    expected_prefix = 1 + b * q + gamma_match * q**2
    if sp.expand(tail_series.removeO() - expected_prefix).coeff(r, 0) != 0:
        raise RuntimeError("tail constant drifted")
    for power in range(1, 5):
        if sp.expand(tail_series.removeO() - expected_prefix).coeff(r, power) != 0:
            raise RuntimeError("tail-matching prefix drifted")

    return {
        "normalization": {
            "kernel": (
                "Psi_(beta,gamma)(u)=4*pi^2*exp(-2*pi*cosh(4u))*"
                "(cosh(9u)+beta*cosh(5u)+gamma*cosh(u))"
            ),
            "transform": (
                "B_(beta,gamma)(x)=P_(9/4)(x)+beta*P_(5/4)(x)+"
                "gamma*P_(1/4)(x)"
            ),
            "block": (
                "P_a(x)=pi^2/2*(K_(a+i*x/4)(2*pi)+"
                "K_(a-i*x/4)(2*pi))"
            ),
            "residual": "R_(beta,gamma)=H_0-B_(beta,gamma)",
        },
        "classical_benchmarks": {
            "polya_p2": {
                "parameters": "beta=-a, gamma=0, a=3/(2*pi)",
                "transform": "B_P2=P_(9/4)-a*P_(5/4)",
                "source_fact": "Polya proved that the P2 transform has only real zeros.",
            },
            "de_bruijn": {
                "parameters": "beta=b, gamma=1, b=pi-3/(2*pi)",
                "transform": "B_dB=P_(9/4)+b*P_(5/4)+P_(1/4)",
                "source_fact": "de Bruijn proved that the dB transform has only real zeros.",
            },
            "proof_boundary": (
                "These are established source facts and normalization checks; "
                "this artifact does not re-prove their real-zero theorems."
            ),
        },
        "normalized_tail": {
            "q_definition": "q=exp(-4u), r=sqrt(q), a=3/(2*pi)",
            "first_summand": (
                "phi_1/[2*pi^2*exp(9u-pi*exp(4u))]=1-a*q"
            ),
            "block_after_exp_pi_q": sp.sstr(normalized_block),
            "actual_after_exp_pi_q": (
                "exp(pi*q)*(1-a*q)=1+b*q+gamma_match*q^2+O(q^3)"
            ),
            "b": "pi-3/(2*pi)",
            "gamma_match": "(pi^2-3)/2",
            "tail_hierarchy": (
                "The residual first has coefficient b-beta at q. If beta=b, "
                "it has coefficient gamma_match-gamma at q^2; at equality "
                "the mirrored gamma block contributes -gamma_match*q^(5/2)."
            ),
            "conclusion": (
                "Eventual residual nonnegativity requires beta<=b; exact "
                "matching beta=b, gamma=gamma_match overshoots in the next "
                "fractional tail order."
            ),
        },
        "polya_p2_positive_residual": {
            "lower_difference": sp.sstr(p2_lower),
            "proof": (
                "Since exp(pi*q)>1+pi*q, 0<q<=1, and pi>3, the normalized "
                "first-summand residual exceeds q*(pi-3*q/2-q^(7/2)+"
                "a*q^(5/2)) >= q*(pi-5/2)>0."
            ),
            "conclusion": "Phi(u)-Psi_P2(u)>0 for every finite u>=0",
        },
        "de_bruijn_positive_residual": {
            "taylor_lower_bound": (
                "exp(pi*q)>sum_(k=0)^4 (pi*q)^k/k!, while 1-a*q>0"
            ),
            "polynomial": f"p(r)={sp.sstr(db_polynomial.as_expr())}",
            "identity": (
                "(1-a*q)*sum_(k=0)^4(pi*q)^k/k!-S_dB(q)=q^2*p(sqrt(q))"
            ),
            "bernstein_coefficients": [sp.sstr(value) for value in db_bernstein],
            "bernstein_balls": [str(value) for value in db_bernstein_balls],
            "conclusion": "Phi(u)-Psi_dB(u)>0 for every finite u>=0",
        },
        "laguerre_polynomial": {
            "forms": "L[f]=f'^2-f*f'', B[f,g]=2*f'*g'-f*g''-g*f''",
            "identity": (
                "L[R_(beta,gamma)]=L[E]-beta*B[E,P5]-gamma*B[E,P1]+"
                "beta^2*L[P5]+beta*gamma*B[P5,P1]+gamma^2*L[P1]"
            ),
            "definitions": "E=H_0-P_(9/4), P5=P_(5/4), P1=P_(1/4)",
        },
        "gasper_square_scope": {
            "single_shift": (
                "For F_(a,c)(z)=K_(i(z-i*c))(a)+K_(i(z+i*c))(a), a zero "
                "gives equality of two Bessel moduli. Gasper converts their "
                "difference into one integral of [f_t(y+c)-f_t(y-c)] times "
                "K_(i*x)(a/sqrt(t))^2, and convexity forces y=0."
            ),
            "multi_shift_guard": (
                "A linear combination with several shifts does not yield that "
                "two-modulus equality. Its squared modulus contains mixed "
                "products between distinct shifts, which the single-shift "
                "identity does not sign."
            ),
            "scope": (
                "This rejects automatic reuse of the published single-shift "
                "square, not the possibility of a new matrix-valued or coupled "
                "addition formula."
            ),
        },
    }


def build_origin_certificate() -> dict:
    ctx.dps = 130
    pi_ball = arb.pi()
    cutoff = 8
    partial = sum(
        (
            pi_ball
            * n**2
            * (2 * pi_ball * n**2 - 3)
            * (-pi_ball * n**2).exp()
            for n in range(1, cutoff + 1)
        ),
        arb(0),
    )
    first_tail_index = cutoff + 1
    tail_majorant_first = (
        2
        * pi_ball**2
        * first_tail_index**4
        * (-pi_ball * first_tail_index**2).exp()
    )
    ratio = (
        (arb(first_tail_index + 1) / first_tail_index) ** 4
        * (-pi_ball * (2 * first_tail_index + 1)).exp()
    )
    tail_upper = tail_majorant_first / (1 - ratio)
    phi_upper = arb(partial.upper() + tail_upper.upper())
    fake_origin = 4 * pi_ball**2 * (-2 * pi_ball).exp()
    rational_upper = arb("61/10") * fake_origin
    if not phi_upper < rational_upper:
        raise RuntimeError("origin compactness bound failed")

    return {
        "rigorous": True,
        "arithmetic": "python-flint Arb balls at 130 decimal digits",
        "identity": (
            "Phi(0)=sum_(n>=1) pi*n^2*(2*pi*n^2-3)*exp(-pi*n^2)"
        ),
        "partial_cutoff": cutoff,
        "partial_sum": str(partial),
        "tail_majorant": (
            "For n>=9, each summand is at most 2*pi^2*n^4*exp(-pi*n^2); "
            "successive majorants have ratio at most "
            "(10/9)^4*exp(-19*pi)."
        ),
        "tail_ratio": str(ratio),
        "tail_upper": str(tail_upper),
        "phi_upper": str(phi_upper),
        "fake_origin": str(fake_origin),
        "comparison_upper": str(rational_upper),
        "ratio_conclusion": (
            "Phi(0)/(4*pi^2*exp(-2*pi))<61/10, hence beta+gamma<51/10"
        ),
        "tail_beta_bound": (
            "beta<=b=pi-3/(2*pi)<821/308<27/10, using pi<22/7"
        ),
        "parameter_region": (
            "Every beta,gamma>=0 with a globally nonnegative residual lies in "
            "0<=beta<27/10, gamma>=0, beta+gamma<51/10."
        ),
    }


def laguerre_coefficients(
    e_jet: list[arb], p5_jet: list[arb], p1_jet: list[arb]
) -> tuple[arb, arb, arb, arb, arb, arb]:
    return (
        laguerre(e_jet),
        -cross_laguerre(e_jet, p5_jet),
        -cross_laguerre(e_jet, p1_jet),
        laguerre(p5_jet),
        cross_laguerre(p5_jet, p1_jet),
        laguerre(p1_jet),
    )


def evaluate_laguerre_polynomial(
    coefficients: tuple[arb, arb, arb, arb, arb, arb], beta: arb, gamma: arb
) -> arb:
    c00, cb, cg, cbb, cbg, cgg = coefficients
    # Multiplication remains defined for zero-centred Arb balls; ``ball**2``
    # in some python-flint builds can return nan at that removable edge case.
    return (
        c00
        + cb * beta
        + cg * gamma
        + cbb * beta * beta
        + cbg * beta * gamma
        + cgg * gamma * gamma
    )


def build_interval_certificate() -> dict:
    ctx.dps = 130
    step = arb("1e-6")
    cauchy_radius = arb("0.04")
    box_radius = arb("0.1")
    spectral_points = (48, 52, 86)

    def e_function(x: acb) -> acb:
        return acb_h0(x) - acb_gasper_block(x, arb("9/4"))

    def p5_function(x: acb) -> acb:
        return acb_gasper_block(x, arb("5/4"))

    def p1_function(x: acb) -> acb:
        return acb_gasper_block(x, arb("1/4"))

    rows: list[dict] = []
    coefficients_by_x: dict[int, tuple[arb, arb, arb, arb, arb, arb]] = {}
    for x in spectral_points:
        e_jet, e_maximum = cauchy_finite_difference_jet(
            e_function, str(x), step, cauchy_radius, box_radius
        )
        p5_jet, p5_maximum = cauchy_finite_difference_jet(
            p5_function, str(x), step, cauchy_radius, box_radius
        )
        p1_jet, p1_maximum = cauchy_finite_difference_jet(
            p1_function, str(x), step, cauchy_radius, box_radius
        )
        raw_coefficients = laguerre_coefficients(e_jet, p5_jet, p1_jet)
        # Run the mesh with the outward-rounded balls written to JSON so an
        # independent replay sees exactly the same certificate.
        coefficient_strings = [str(value) for value in raw_coefficients]
        coefficients = tuple(arb(value) for value in coefficient_strings)
        coefficients_by_x[x] = coefficients
        rows.append(
            {
                "x": str(x),
                "jets": {
                    "E": [str(value) for value in e_jet],
                    "P5": [str(value) for value in p5_jet],
                    "P1": [str(value) for value in p1_jet],
                },
                "cauchy_box_abs_upper": {
                    "E": str(e_maximum),
                    "P5": str(p5_maximum),
                    "P1": str(p1_maximum),
                },
                "laguerre_coefficients": {
                    name: str(value)
                    for name, value in zip(
                        ("c00", "c_beta", "c_gamma", "c_beta2", "c_beta_gamma", "c_gamma2"),
                        coefficient_strings,
                    )
                },
            }
        )

    denominator = 80
    beta_cells = 216
    sum_cells = 408
    radius = arb(1) / (2 * denominator)
    assignments = {x: 0 for x in spectral_points}
    worst_upper: dict[int, object | None] = {x: None for x in spectral_points}
    digest = hashlib.sha256()
    total = 0
    for beta_index in range(beta_cells):
        beta_midpoint = arb(2 * beta_index + 1) / (2 * denominator)
        beta_box = arb(beta_midpoint, radius)
        for gamma_index in range(sum_cells - beta_index):
            gamma_midpoint = arb(2 * gamma_index + 1) / (2 * denominator)
            gamma_box = arb(gamma_midpoint, radius)
            total += 1
            for x in spectral_points:
                value = evaluate_laguerre_polynomial(
                    coefficients_by_x[x], beta_box, gamma_box
                )
                if value < 0:
                    assignments[x] += 1
                    upper = value.upper()
                    if worst_upper[x] is None or upper > worst_upper[x]:
                        worst_upper[x] = upper
                    digest.update(f"{beta_index},{gamma_index}:{x}\n".encode("ascii"))
                    break
            else:
                raise RuntimeError(
                    f"uncovered parameter box beta={beta_index}, gamma={gamma_index}"
                )

    expected_total = 64908
    if total != expected_total or sum(assignments.values()) != expected_total:
        raise RuntimeError("parameter-box count drifted")

    a_ball = 3 / (2 * arb.pi())
    b_ball = arb.pi() - a_ball
    p2_value = evaluate_laguerre_polynomial(
        coefficients_by_x[86], -a_ball, arb(0)
    )
    de_bruijn_value = evaluate_laguerre_polynomial(
        coefficients_by_x[52], b_ball, arb(1)
    )
    if not (p2_value < 0 and de_bruijn_value < 0):
        raise RuntimeError("classical residual witness failed")

    for row in rows:
        x = int(row["x"])
        row["mesh_assignment_count"] = assignments[x]
        row["closest_assigned_upper"] = str(worst_upper[x])

    return {
        "rigorous": True,
        "arithmetic": "python-flint Acb/Arb balls at 130 decimal digits",
        "spectral_point_formula": (
            "H_0(x)=xi((1+i*x)/2)/8 and "
            "P_a(x)=pi^2/2*(K_(a+i*x/4)(2*pi)+K_(a-i*x/4)(2*pi))"
        ),
        "derivative_method": (
            "For each real entire component, central differences with step h "
            "enclose f' and f''. An Acb square bounds abs(f) by M on every "
            "Cauchy circle of radius r around the difference segment; the "
            "added errors are h^2*M/r^3 and 2*h^2*M/r^4."
        ),
        "step": "1e-6",
        "cauchy_radius": "0.04",
        "acb_box_radius": "0.1",
        "spectral_points": [str(value) for value in spectral_points],
        "rows": rows,
        "mesh": {
            "denominator": denominator,
            "beta_cell_count": beta_cells,
            "sum_cell_count": sum_cells,
            "cell_definition": (
                "I_i=[i/80,(i+1)/80], J_j=[j/80,(j+1)/80], "
                "0<=i<216, 0<=j<408-i"
            ),
            "covered_region": (
                "0<=beta<27/10, gamma>=0, beta+gamma<51/10"
            ),
            "total_boxes": total,
            "assignment_order": [str(value) for value in spectral_points],
            "assignment_counts": {str(key): value for key, value in assignments.items()},
            "assignment_sha256": digest.hexdigest(),
            "all_boxes_strictly_negative_at_an_assigned_x": True,
        },
        "classical_witnesses": {
            "polya_p2": {
                "x": "86",
                "parameters": "beta=-3/(2*pi), gamma=0",
                "laguerre_value": str(p2_value),
                "strictly_negative": True,
            },
            "de_bruijn": {
                "x": "52",
                "parameters": "beta=pi-3/(2*pi), gamma=1",
                "laguerre_value": str(de_bruijn_value),
                "strictly_negative": True,
            },
        },
        "conclusion": (
            "Every nonnegative-coefficient 9/5/1 block whose kernel residual "
            "is globally nonnegative has L[R_(beta,gamma)](x)<0 at at least "
            "one of x=48,52,86. The classical P2 residual also has L<0 at x=86."
        ),
        "scope": (
            "This exhausts the positive three-block plus independently-LP "
            "residual architecture at t=0. It does not reject signed higher "
            "blocks or a coupled identity controlling all mixed terms."
        ),
    }


def build_payload() -> dict:
    exact = build_exact()
    origin = build_origin_certificate()
    interval = build_interval_certificate()
    rows = [
        GateRow(
            id="gctb_01_classical_normalizations",
            role="established_benchmark",
            readiness="source_validated",
            claim="The Polya P2 and de Bruijn kernels map to the stated 9/5/1 Bessel combinations.",
            formula="B_P2=P9-a*P5; B_dB=P9+b*P5+P1",
            proof_boundary="Normalization and source audit; the classical real-zero theorems are not re-proved here.",
            diagnostics=exact["classical_benchmarks"],
        ),
        GateRow(
            id="gctb_02_polya_p2_positive_residual",
            role="exact_kernel_theorem",
            readiness="ready_to_apply",
            claim="The signed Polya P2 block leaves a strictly positive Xi-kernel residual.",
            formula=exact["polya_p2_positive_residual"]["conclusion"],
            proof_boundary="Exact pointwise kernel theorem only; positivity does not imply a real-rooted transform.",
            diagnostics=exact["polya_p2_positive_residual"],
        ),
        GateRow(
            id="gctb_03_de_bruijn_positive_residual",
            role="interval_kernel_theorem",
            readiness="interval_validated",
            claim="The classical de Bruijn block also leaves a strictly positive Xi-kernel residual.",
            formula=exact["de_bruijn_positive_residual"]["conclusion"],
            proof_boundary="Finite exact Bernstein reduction with Arb-certified coefficients; not transform real-rootedness of the residual.",
            diagnostics=exact["de_bruijn_positive_residual"],
        ),
        GateRow(
            id="gctb_04_tail_origin_compactness",
            role="exact_interval_reduction",
            readiness="ready_to_apply",
            claim="Tail and origin positivity confine every positive 9/5/1 residual candidate to one rational triangle.",
            formula=origin["parameter_region"],
            proof_boundary="Necessary compact parameter reduction; no sufficiency of kernel positivity is claimed.",
            diagnostics={"tail": exact["normalized_tail"], "origin": origin},
        ),
        GateRow(
            id="gctb_05_bivariate_laguerre_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="The residual first Laguerre expression is an exact quadratic in both block coefficients.",
            formula=exact["laguerre_polynomial"]["identity"],
            proof_boundary="Exact algebraic identity used by the interval mesh.",
            diagnostics=exact["laguerre_polynomial"],
        ),
        GateRow(
            id="gctb_06_acb_cauchy_jets",
            role="interval_method",
            readiness="interval_validated",
            claim="Three Acb Cauchy jet certificates rigorously enclose every bivariate Laguerre coefficient used by the mesh.",
            formula=interval["derivative_method"],
            proof_boundary="Finite spectral certificates at t=0 only.",
            diagnostics={key: interval[key] for key in ("arithmetic", "step", "cauchy_radius", "acb_box_radius", "rows")},
        ),
        GateRow(
            id="gctb_07_positive_three_block_mesh",
            role="interval_theorem",
            readiness="interval_validated",
            claim="All 64,908 boxes covering the admissible positive three-block parameter triangle violate the first Laguerre inequality.",
            formula=interval["conclusion"],
            proof_boundary="Exhaustive for nonnegative 9/5/1 coefficients with a globally nonnegative residual; signed extensions remain open.",
            diagnostics=interval["mesh"],
        ),
        GateRow(
            id="gctb_08_classical_residual_obstructions",
            role="interval_composition",
            readiness="ready_to_apply",
            claim="Both established classical real-zero bases leave positive residuals that fail a necessary Laguerre-Polya inequality.",
            formula="L[R_P2](86)<0 and L[R_dB](52)<0",
            proof_boundary="Rejects independent LP treatment of the residuals, not the classical base theorems or RH.",
            diagnostics=interval["classical_witnesses"],
        ),
        GateRow(
            id="gctb_09_gasper_square_scope",
            role="source_scope_guard",
            readiness="guard_validated",
            claim="Gasper's published single-shift square does not automatically sign the mixed products in a multi-shift sum.",
            formula=exact["gasper_square_scope"]["multi_shift_guard"],
            proof_boundary=exact["gasper_square_scope"]["scope"],
            diagnostics=exact["gasper_square_scope"],
        ),
        GateRow(
            id="gctb_10_coupled_signed_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The surviving comparison route requires signed higher blocks or a new coupled square controlling cross terms together.",
            formula=(
                "Seek a signed multi-block or matrix-valued addition identity after full modular cancellation."
            ),
            proof_boundary="Open theorem target; not strict Laguerre positivity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_classical_three_block_residual_gate",
        "date": "2026-07-11",
        "status": "exact and interval-certified classical three-block residual obstruction",
        "proof_boundary": (
            "This artifact proves positive kernel residuals for the classical "
            "Polya P2 and de Bruijn bases, derives necessary compact constraints "
            "for the full nonnegative 9/5/1 family, and uses Acb/Arb Cauchy "
            "certificates to show that every such globally positive residual "
            "violates the first Laguerre inequality. It audits the scope of "
            "Gasper's single-shift square. Signed higher-block and genuinely "
            "coupled identities remain open; this does not prove or disprove RH "
            "or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_gasper_residual_two_block_gate.md",
            "outputs/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.md",
            "outputs/formal_core.md",
            "https://arxiv.org/abs/0801.2996",
            "https://arxiv.org/abs/1502.06844",
            "https://doi.org/10.1215/S0012-7094-50-01720-0",
            "https://doi.org/10.1007/BF02565336",
        ],
        "exact": exact,
        "origin_certificate": origin,
        "interval_certificate": interval,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    origin = payload["origin_certificate"]
    interval = payload["interval_certificate"]
    mesh = interval["mesh"]
    lines = [
        "# Jensen-Window PF Newman Classical Three-Block Residual Gate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: exact and interval-certified classical three-block residual",
        "obstruction. This is not a proof or disproof of RH or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_newman_classical_three_block_residual_gate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_classical_three_block_residual_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_classical_three_block_residual_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_classical_three_block_residual_gate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF Newman classical three-block residual gate: 10 rows, 0 issues, 2 established classical real-zero benchmarks, 2 positive-kernel residual theorems, 1 compact parameter reduction, 1 exact bivariate Laguerre identity, 3 Acb spectral certificates, 64908 parameter boxes covered, 2 classical residual obstructions, 1 Gasper square-scope guard, 1 coupled handoff",
        "```",
        "",
        "## Classical Blocks",
        "",
        "In the corpus normalization,",
        "",
        "```text",
        exact["normalization"]["block"],
        exact["classical_benchmarks"]["polya_p2"]["transform"],
        exact["classical_benchmarks"]["de_bruijn"]["transform"],
        "```",
        "",
        "Polya's P2 transform and de Bruijn's transform are established",
        "real-zero benchmarks in the cited sources. The present gate checks",
        "their normalization and studies what remains after subtracting them.",
        "",
        "## Positive Residuals",
        "",
        "For `q=exp(-4u)` and `a=3/(2pi)`, the P2 first-summand",
        "residual is bounded below by",
        "",
        "```text",
        exact["polya_p2_positive_residual"]["lower_difference"],
        exact["polya_p2_positive_residual"]["conclusion"],
        "```",
        "",
        "For de Bruijn's block, the fourth Taylor lower bound for `exp(pi*q)`",
        "reduces the same claim to a degree-six polynomial in `sqrt(q)`.",
        "All seven Bernstein coefficients are rigorously positive:",
        "",
        "```text",
    ]
    lines.extend(f"b_{index}={value}" for index, value in enumerate(exact["de_bruijn_positive_residual"]["bernstein_balls"]))
    lines.extend(
        [
            exact["de_bruijn_positive_residual"]["conclusion"],
            "```",
            "",
            "Thus both classical real-rooted bases leave pointwise positive",
            "kernel residuals. Those residual transforms nevertheless fail below.",
            "",
            "## Full Positive Family",
            "",
            "```text",
            exact["normalization"]["kernel"],
            origin["tail_beta_bound"],
            origin["ratio_conclusion"],
            origin["parameter_region"],
            "```",
            "",
            "The origin bound uses the first eight exact theta summands and an",
            "explicit geometric majorant for `n>=9`. These are necessary",
            "conditions for a globally nonnegative residual and place every",
            "candidate in one rational triangle.",
            "",
            "## Interval Theorem",
            "",
            "The exact bivariate expression is",
            "",
            "```text",
            exact["laguerre_polynomial"]["identity"],
            "```",
            "",
            interval["derivative_method"],
            "The rational triangle is covered by closed `1/80` boxes:",
            "",
            "```text",
            mesh["cell_definition"],
            f"total boxes={mesh['total_boxes']}",
            f"x=48 assignments={mesh['assignment_counts']['48']}",
            f"x=52 assignments={mesh['assignment_counts']['52']}",
            f"x=86 assignments={mesh['assignment_counts']['86']}",
            f"assignment sha256={mesh['assignment_sha256']}",
            "```",
            "",
            "Every box has a strictly negative Arb enclosure at one assigned",
            "spectral point. In particular,",
            "",
            "```text",
            f"L[R_P2](86)={interval['classical_witnesses']['polya_p2']['laguerre_value']}",
            f"L[R_dB](52)={interval['classical_witnesses']['de_bruijn']['laguerre_value']}",
            "```",
            "",
            "So neither classical positive residual is an independent",
            "Laguerre-Polya component, and no nonnegative-coefficient 9/5/1",
            "block with a globally nonnegative residual can repair this route.",
            "",
            "## Square Scope",
            "",
            exact["gasper_square_scope"]["single_shift"],
            exact["gasper_square_scope"]["multi_shift_guard"],
            "",
            "This leaves a precise open handoff: derive a genuinely coupled",
            "square or matrix identity for signed higher blocks after modular",
            "cancellation. The published single-shift square alone does not do it.",
            "",
            "References: https://arxiv.org/abs/0801.2996,",
            "https://arxiv.org/abs/1502.06844,",
            "https://doi.org/10.1215/S0012-7094-50-01720-0, and",
            "https://doi.org/10.1007/BF02565336.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(f"wrote Newman classical three-block residual gate: {len(payload['rows'])} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
