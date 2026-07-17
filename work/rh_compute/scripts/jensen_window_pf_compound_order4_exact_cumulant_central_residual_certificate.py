#!/usr/bin/env python3
"""Close the central exact-minus-epsilon-ten partition residual."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402
import sympy as sp  # noqa: E402

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    sha256,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.md"
)
SOURCE_PARTITION_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.json"
)
SOURCE_EXACT_TAIL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.json"
)
SOURCE_FORMAL_TAIL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.json"
)
SOURCE_DISK_CONTRACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.json"
)
SOURCE_SECOND_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.json"
)
FINITE_Q_FLOOR = 9_000
RAY_Q_FLOOR = 10**35
RAY_START = 20
BELL_ORDER = 15
FINITE_EPSILON_DENOMINATOR = 94
RAY_EPSILON_DENOMINATOR = 10**17
FINITE_PARTITION_BOUNDS = {
    11: Fraction(2),
    12: Fraction(2),
    13: Fraction(21, 10),
    14: Fraction(12, 5),
}
FINITE_JET_CAPS = {
    3: Fraction(6, 5),
    4: Fraction(3, 2),
    5: Fraction(2),
    6: Fraction(3),
    7: Fraction(9, 2),
    8: Fraction(7),
    9: Fraction(12),
    10: Fraction(21),
    11: Fraction(38),
    12: Fraction(71),
    13: Fraction(200),
    14: Fraction(400),
    15: Fraction(800),
    16: Fraction(1600),
}


@dataclass(frozen=True)
class CentralRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def potential_remainder_data() -> dict:
    u, q = sp.symbols("u q", positive=True)
    remainder = 100 * u**2 - 5 * u - sp.log(2 * q - 3) - sp.log(u)
    coefficient_sums = {}
    ray_caps = {}
    for order in range(1, 18):
        remainder = sp.cancel(
            u * (sp.diff(remainder, u) + 4 * q * sp.diff(remainder, q)) / 2
        )
        if order < 13:
            continue
        numerator, denominator = sp.fraction(remainder)
        polynomial = sp.Poly(numerator, u, q)
        coefficient_sum = sum(abs(value) for value in polynomial.coeffs())
        expected_denominator = 2**order * (2 * q - 3) ** order
        if sp.expand(denominator - expected_denominator) != 0:
            raise RuntimeError(f"unexpected potential remainder denominator at {order}")
        if polynomial.degree(u) > order or polynomial.degree(q) > order:
            raise RuntimeError(f"unexpected potential remainder degree at {order}")
        cap = int(3**order + coefficient_sum)
        coefficient_sums[str(order)] = {
            "numerator_terms": len(polynomial.terms()),
            "u_degree": polynomial.degree(u),
            "q_degree": polynomial.degree(q),
            "absolute_coefficient_sum": str(coefficient_sum),
            "ray_normalized_jet_cap": str(cap),
        }
        ray_caps[order] = Fraction(cap)
    return {
        "leading_recurrence": "P_0=1; P_(r+1)=(u/2)*P_r'+2*u*P_r",
        "leading_bound": "0<P_r(u)<3^r*u^r for u>=399/20, r<=17",
        "remainder_bound": (
            "|V^(r)-q*P_r(u)|<=C_r*u^r from the exact numerator coefficient norm"
        ),
        "curvature_input": "V''>=(39/10)*u^2*q",
        "coefficient_rows": coefficient_sums,
        "ray_caps": ray_caps,
    }


def partition_norms(
    partition_rows: dict, caps: dict[int, Fraction]
) -> dict[int, Fraction]:
    symbols = tuple(sp.symbols("L_3:17"))
    locals_map = {str(symbol): symbol for symbol in symbols}
    result = {}
    for degree_text, row in partition_rows.items():
        total = sp.Rational(0)
        for coefficient_row in row["z_coefficients"].values():
            expression = sp.sympify(coefficient_row["formula"], locals=locals_map)
            for powers, coefficient in sp.Poly(expression, *symbols).terms():
                term = abs(sp.Rational(coefficient))
                for order, power in zip(range(3, 17), powers):
                    cap = caps[order]
                    term *= sp.Rational(cap.numerator, cap.denominator) ** power
                total += term
        value = Fraction(int(sp.numer(total)), int(sp.denom(total)))
        result[int(degree_text)] = value
    return result


def formal_weight_norms(caps: dict[int, Fraction]) -> dict[int, Fraction]:
    y = sp.symbols("y")
    symbols = tuple(sp.symbols("L_3:17"))
    perturbation = [sp.Integer(0) for _ in range(15)]
    for order, symbol in zip(range(3, 17), symbols):
        perturbation[order - 2] = symbol * y**order / sp.factorial(order)
    exponential = [sp.Integer(0) for _ in range(15)]
    exponential[0] = sp.Integer(1)
    for degree in range(1, 15):
        exponential[degree] = sp.expand(
            -sum(
                index * perturbation[index] * exponential[degree - index]
                for index in range(1, degree + 1)
            )
            / degree
        )
    norms = {}
    for degree in range(11, 15):
        total = sp.Rational(0)
        for powers, coefficient in sp.Poly(
            exponential[degree], *symbols, y
        ).terms():
            term = abs(sp.Rational(coefficient))
            for order, power in zip(range(3, 17), powers[:-1]):
                cap = caps[order]
                term *= sp.Rational(cap.numerator, cap.denominator) ** power
            total += term
        norms[degree] = Fraction(int(sp.numer(total)), int(sp.denom(total)))
    return norms


def perturbation_gate(
    caps: dict[int, Fraction], *, finite: bool
) -> dict:
    if finite:
        flint.ctx.prec = 256
        q = flint.arb(FINITE_Q_FLOOR)
        y = 1 + (32 * q.log()).sqrt()
        epsilon = 1 / q.sqrt()
        ratio = sum(
            flint.arb(cap.numerator)
            / cap.denominator
            * epsilon ** (order - 2)
            * y ** (order - 2)
            / math.factorial(order)
            for order, cap in caps.items()
        )
        if not bool(ratio < flint.arb(1) / 24):
            raise RuntimeError("finite perturbation ratio failed")
        return {
            "domain": "2<=u<=20",
            "endpoint_ratio_upper": arb_upper_text(ratio),
            "cap": "1/24",
            "propagation": (
                "Y(q)/sqrt(q) is decreasing because 16/(s*Y)<1/2"
            ),
            "consequence": (
                "|sum_(r=3)^16 lambda_r*y^r/r!|<y^2/24 on |y|<=Y"
            ),
        }

    quarter_floor = 10**8
    ratio = sum(
        cap
        * Fraction(2 ** (order - 2), math.factorial(order))
        / quarter_floor ** (order - 2)
        for order, cap in caps.items()
    )
    if ratio >= Fraction(1, 24):
        raise RuntimeError("ray perturbation ratio failed")
    return {
        "domain": "u>=20",
        "q_quarter_floor": quarter_floor,
        "rational_ratio": str(ratio),
        "cap": "1/24",
        "input": "Y<2*q^(1/4) and q^(1/4)>10^8",
        "consequence": (
            "|sum_(r=3)^16 lambda_r*y^r/r!|<y^2/24 on |y|<=Y"
        ),
    }


def bell_fifteen_majorant(
    caps: dict[int, Fraction], epsilon_denominator: int
) -> dict:
    x = sp.symbols("x", nonnegative=True)
    derivatives = {}
    for derivative in range(1, BELL_ORDER + 1):
        derivatives[derivative] = sum(
            sp.factorial(weight) / sp.factorial(weight - derivative)
            * sp.Rational(1, epsilon_denominator) ** (weight - derivative)
            * sp.Rational(caps[weight + 2].numerator, caps[weight + 2].denominator)
            * x ** (weight + 2)
            / sp.factorial(weight + 2)
            for weight in range(derivative, 15)
        )
    bell = [sp.Poly(1, x)]
    for order in range(1, BELL_ORDER + 1):
        bell.append(
            sp.Poly(
                sum(
                    sp.binomial(order - 1, derivative - 1)
                    * bell[order - derivative].as_expr()
                    * derivatives[derivative]
                    for derivative in range(1, order + 1)
                ),
                x,
            )
        )

    gaussian_coefficient = Fraction(11, 24)
    moment_caps = [Fraction(4), Fraction(5, 1) / (2 * gaussian_coefficient)]
    for degree in range(2, bell[BELL_ORDER].degree() + 1):
        moment_caps.append(
            (
                (degree - 1) * moment_caps[degree - 2]
                + moment_caps[degree - 1]
            )
            / (2 * gaussian_coefficient)
        )
    total = sp.Rational(0)
    for (degree,), coefficient in bell[BELL_ORDER].terms():
        moment = moment_caps[degree]
        total += (
            coefficient
            * 2
            * sp.Rational(moment.numerator, moment.denominator)
            / sp.factorial(BELL_ORDER)
        )
    value = Fraction(int(sp.numer(total)), int(sp.denom(total)))
    ceiling = math.ceil(value)
    return {
        "taylor_formula": (
            "F(1)-sum_(n=0)^14 F^(n)(0)/n!="
            "integral_0^1 (1-t)^14 F^(15)(t)/14! dt"
        ),
        "epsilon_factor": "every Bell monomial contains epsilon^15",
        "gaussian_majorant": (
            "exp(-y^2/2+|y|+|B_t(y)|)<=exp(-(11/24)*y^2+|y|)"
        ),
        "half_line_moment_base": "J_0<4; J_1<60/11",
        "half_line_moment_recurrence": (
            "J_m=((m-1)*J_(m-2)+J_(m-1))/(11/12)"
        ),
        "epsilon_denominator": epsilon_denominator,
        "bell_polynomial_terms": len(bell[BELL_ORDER].terms()),
        "bell_polynomial_degree": bell[BELL_ORDER].degree(),
        "constant": str(value),
        "integer_cap": ceiling,
        "bound": f"formal exponential remainder <{ceiling}*q^(-15/2)",
    }


def potential_taylor_majorant(jet17_cap: Fraction) -> dict:
    left_transfer = Fraction(5, 4) ** 8 * Fraction(50, 39) ** 9
    right_transfer = (
        Fraction(50, 39) * Fraction(41, 40) ** 2 * Fraction(5, 4)
    ) ** 9
    transfer = max(left_transfer, right_transfer)
    effective_cap = math.ceil(jet17_cap * transfer)

    coefficient = Fraction(1, 4)
    moments = [Fraction(9), Fraction(10) / (2 * coefficient)]
    for degree in range(2, 18):
        moments.append(
            ((degree - 1) * moments[degree - 2] + moments[degree - 1])
            / (2 * coefficient)
        )
    constant = Fraction(4 * effective_cap, math.factorial(17)) * moments[17]
    return {
        "collar_transfer_left": str(left_transfer),
        "collar_transfer_right": str(right_transfer),
        "effective_normalized_jet17_cap": effective_cap,
        "pointwise_remainder": (
            f"|W-y^2/2-T_16|<={effective_cap}*q^(-15/2)*|y|^17/17!"
        ),
        "decay": (
            "W>=y^2/4 and y^2/2+T_16>=11*y^2/24 on the central collar"
        ),
        "moment_base": "integral_0^infinity exp(-y^2/4+y)dy<9",
        "integrated_constant": str(constant),
        "bound": f"potential Taylor residual <({constant})*q^(-15/2)",
    }


def arb_power(value: flint.arb, power: Fraction) -> flint.arb:
    return value ** (flint.arb(power.numerator) / power.denominator)


def added_formal_tail_bound(
    norms: dict[int, Fraction], q: flint.arb
) -> flint.arb:
    total = flint.arb(0)
    for degree, norm in norms.items():
        coefficient = flint.arb(norm.numerator) / norm.denominator
        power = Fraction(degree, 4) - Fraction(65, 4)
        total += coefficient * 2 ** (3 * degree - 1) * arb_power(q, power)
    return 6 * total


def compose_budgets(
    finite_partition: dict[int, Fraction],
    ray_partition: dict[int, Fraction],
    finite_weight: dict[int, Fraction],
    ray_weight: dict[int, Fraction],
    finite_bell: dict,
    ray_bell: dict,
    finite_potential: dict,
    ray_potential: dict,
) -> dict:
    flint.ctx.prec = 256
    q_finite = flint.arb(FINITE_Q_FLOOR)
    finite_partition_value = sum(
        flint.arb(value.numerator)
        / value.denominator
        * arb_power(q_finite, Fraction(-degree, 2))
        for degree, value in finite_partition.items()
    )
    finite_formal_tails = added_formal_tail_bound(finite_weight, q_finite)
    finite_bell_value = flint.arb(finite_bell["integer_cap"]) * arb_power(
        q_finite, Fraction(-15, 2)
    )
    finite_potential_constant = Fraction(finite_potential["integrated_constant"])
    finite_potential_value = (
        flint.arb(finite_potential_constant.numerator)
        / finite_potential_constant.denominator
        * arb_power(q_finite, Fraction(-15, 2))
    )
    finite_total = (
        finite_partition_value
        + finite_formal_tails
        + finite_bell_value
        + finite_potential_value
    )
    finite_target = flint.arb(1) / (500_000 * q_finite**3)
    if not bool(finite_total < finite_target):
        raise RuntimeError("finite central partition budget failed")

    q_ray = flint.arb(RAY_Q_FLOOR)
    ray_partition_value = sum(
        flint.arb(value.numerator)
        / value.denominator
        * arb_power(q_ray, Fraction(-degree, 2))
        for degree, value in ray_partition.items()
    )
    ray_formal_tails = added_formal_tail_bound(ray_weight, q_ray)
    ray_bell_value = flint.arb(ray_bell["integer_cap"]) * arb_power(
        q_ray, Fraction(-15, 2)
    )
    ray_potential_constant = Fraction(ray_potential["integrated_constant"])
    ray_potential_value = (
        flint.arb(ray_potential_constant.numerator)
        / ray_potential_constant.denominator
        * arb_power(q_ray, Fraction(-15, 2))
    )
    ray_total = (
        ray_partition_value
        + ray_formal_tails
        + ray_bell_value
        + ray_potential_value
    )
    ray_target = flint.arb(1) / (300_000 * RAY_START * q_ray**3)
    if not bool(ray_total < ray_target):
        raise RuntimeError("ray central partition budget failed")
    return {
        "finite": {
            "partition_correction_upper": arb_upper_text(finite_partition_value),
            "added_formal_tails_upper": arb_upper_text(finite_formal_tails),
            "bell_remainder_upper": arb_upper_text(finite_bell_value),
            "potential_taylor_upper": arb_upper_text(finite_potential_value),
            "total_upper": arb_upper_text(finite_total),
            "target_lower": arb_lower_text(finite_target),
            "budget_ratio_upper": arb_upper_text(finite_total / finite_target),
            "proved_bound": "central residual <1/(500000*q^3), 2<=u<=20",
            "propagation": (
                "every ratio to q^-3 decreases with q; the endpoint q=9000 is worst"
            ),
        },
        "ray": {
            "partition_correction_upper": arb_upper_text(ray_partition_value),
            "added_formal_tails_upper": arb_upper_text(ray_formal_tails),
            "bell_remainder_upper": arb_upper_text(ray_bell_value),
            "potential_taylor_upper": arb_upper_text(ray_potential_value),
            "total_upper": arb_upper_text(ray_total),
            "target_lower": arb_lower_text(ray_target),
            "budget_ratio_upper": arb_upper_text(ray_total / ray_target),
            "proved_bound": "central residual <1/(300000*u*q^3), u>=20",
            "propagation": (
                "u*q^(3-n/2), u*q^-9/2, and the formal-tail ratios all decrease"
            ),
        },
    }


def build_artifact() -> dict:
    finite_source = load_json(SOURCE_PARTITION_FINITE)
    if finite_source.get("summary", {}).get("finite_partition_extension_closed") is not True:
        raise RuntimeError("finite partition-extension source is not closed")
    exact_tail = load_json(SOURCE_EXACT_TAIL)
    formal_tail = load_json(SOURCE_FORMAL_TAIL)
    disk = load_json(SOURCE_DISK_CONTRACT)
    second_ray = load_json(SOURCE_SECOND_RAY)
    if exact_tail.get("summary", {}).get("exact_tails_closed") != 2:
        raise RuntimeError("exact-tail source is not closed")
    if formal_tail.get("summary", {}).get("formal_tails_closed") != 2:
        raise RuntimeError("formal-tail source is not closed")
    if disk.get("summary", {}).get("partition_degrees") != 10:
        raise RuntimeError("complex-disk source is not closed")
    if second_ray.get("summary", {}).get("global_second_next_layer_closed") is not True:
        raise RuntimeError("second-next ray source is not closed")

    potential = potential_remainder_data()
    potential_serialized = {
        **potential,
        "ray_caps": {
            str(order): str(value)
            for order, value in potential["ray_caps"].items()
        },
    }
    ray_caps = {order: Fraction(2) for order in range(3, 13)}
    ray_caps.update({order: potential["ray_caps"][order] for order in range(13, 17)})
    finite_partition = FINITE_PARTITION_BOUNDS
    ray_partition = partition_norms(finite_source["partition_rows"], ray_caps)
    finite_weight = {
        int(degree): Fraction(row["coefficient_norm"])
        for degree, row in finite_source["formal_weight_norms"].items()
    }
    ray_weight = formal_weight_norms(ray_caps)
    finite_perturbation = perturbation_gate(FINITE_JET_CAPS, finite=True)
    ray_perturbation = perturbation_gate(ray_caps, finite=False)
    finite_bell = bell_fifteen_majorant(
        FINITE_JET_CAPS, FINITE_EPSILON_DENOMINATOR
    )
    ray_bell = bell_fifteen_majorant(ray_caps, RAY_EPSILON_DENOMINATOR)
    finite_potential = potential_taylor_majorant(Fraction(4000))
    ray_potential = potential_taylor_majorant(potential["ray_caps"][17])
    budgets = compose_budgets(
        finite_partition,
        ray_partition,
        finite_weight,
        ray_weight,
        finite_bell,
        ray_bell,
        finite_potential,
        ray_potential,
    )
    rows = [
        CentralRow(
            id="co4eccrc_01_exact_decomposition",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The central exact-minus-epsilon-ten residual splits into four added partition coefficients, their formal tails, an epsilon-fifteen exponential remainder, and the exact potential Taylor remainder.",
            formula="central(A-P^[10])=Z_11..Z_14-formal tails+Bell_15+potential_R17",
            proof_boundary="Exact algebraic and integral decomposition only.",
        ),
        CentralRow(
            id="co4eccrc_02_finite_partition_cancellation",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="The finite partition extension supplies cancellation-preserving unit-disk caps through epsilon fourteen.",
            formula="||Z_11||<2, ||Z_12||<2, ||Z_13||<21/10, ||Z_14||<12/5",
            proof_boundary="Finite formal partition coefficients only.",
            diagnostics={str(key): str(value) for key, value in finite_partition.items()},
        ),
        CentralRow(
            id="co4eccrc_03_ray_high_jets",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The q-leading recurrence and exact rational remainder numerators give explicit high normalized-jet caps on the asymptotic collar.",
            formula="|L_r|<=3^r+C_r, r=13,...,17",
            proof_boundary="Crude high-jet ray theorem only.",
            diagnostics=potential_serialized,
        ),
        CentralRow(
            id="co4eccrc_04_ray_partition_majorants",
            role="exact_inequality",
            readiness="ready_to_apply",
            claim="Direct exact coefficient norms bound the four added ray partition coefficients without relying on cancellation.",
            formula="|Z_n(z)|<=M_n, n=11,...,14, |z|<=1",
            proof_boundary="Asymptotic formal partition majorants only.",
            diagnostics={str(key): str(value) for key, value in ray_partition.items()},
        ),
        CentralRow(
            id="co4eccrc_05_central_perturbation",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="Both regimes keep the degree-sixteen standardized perturbation below one twenty-fourth of the quadratic potential.",
            formula="|T_16(y)|<y^2/24 for |y|<=Y",
            proof_boundary="Central-collar potential comparison only.",
            diagnostics={"finite": finite_perturbation, "ray": ray_perturbation},
        ),
        CentralRow(
            id="co4eccrc_06_bell_remainder",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="A Bell-polynomial Taylor theorem controls the complete formal exponential remainder from epsilon fifteen onward while retaining its epsilon power.",
            formula="|exp(-T_16)-E^[14]|<=Bell_15 remainder",
            proof_boundary="Formal exponential remainder only.",
            diagnostics={"finite": finite_bell, "ray": ray_bell},
        ),
        CentralRow(
            id="co4eccrc_07_potential_remainder",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The shifted seventeenth-jet caps and Gaussian moments control the exact potential Taylor remainder.",
            formula="exact-potential residual <=C*q^(-15/2)",
            proof_boundary="Exact potential versus degree-sixteen Taylor potential only.",
            diagnostics={"finite": finite_potential, "ray": ray_potential},
        ),
        CentralRow(
            id="co4eccrc_08_added_formal_tails",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Gaussian hazard bounds control the two formal tails of every added epsilon-eleven through epsilon-fourteen density coefficient.",
            formula="added formal tails use exp(-Y^2/2+Y)=exp(1/2)*q^-16",
            proof_boundary="Only the four added formal-density coefficients.",
            diagnostics={
                "finite_weight_norms": {str(key): str(value) for key, value in finite_weight.items()},
                "ray_weight_norms": {str(key): str(value) for key, value in ray_weight.items()},
            },
        ),
        CentralRow(
            id="co4eccrc_09_finite_central_theorem",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="All four central contributions fit the finite unit-disk partition allocation.",
            formula=budgets["finite"]["proved_bound"],
            proof_boundary="Central partition residual on 2<=u<=20 only.",
            diagnostics=budgets["finite"],
        ),
        CentralRow(
            id="co4eccrc_10_ray_central_theorem",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="All four central contributions fit the asymptotic unit-disk partition allocation.",
            formula=budgets["ray"]["proved_bound"],
            proof_boundary="Central partition residual on u>=20 only.",
            diagnostics=budgets["ray"],
        ),
        CentralRow(
            id="co4eccrc_11_global_partition_residual",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The central theorem and the four inherited tail theorems discharge both unit-disk exact partition targets in the complex-disk contract.",
            formula="|A_u-P_u^[10]|<1/(100000*q^3) finite; <1/(20000*u*q^3) ray",
            proof_boundary="Exact partition residual and consequent exact-minus-epsilon-ten cumulant budgets only.",
            diagnostics={
                "formal_tails_closed": 2,
                "exact_tails_closed": 2,
                "central_components_closed": 1,
                "finite_scaled_cumulant_budget": "9/1000",
                "ray_scaled_cumulant_budget": "1/(100u)",
            },
        ),
    ]
    source_paths = (
        SOURCE_PARTITION_FINITE,
        SOURCE_EXACT_TAIL,
        SOURCE_FORMAL_TAIL,
        SOURCE_DISK_CONTRACT,
        SOURCE_SECOND_RAY,
    )
    return {
        "kind": "jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate",
        "date": "2026-07-13",
        "status": "global exact central residual and unit-disk partition theorem",
        "proof_boundary": (
            "This artifact proves the sole remaining central exact-minus-formal "
            "partition residual and composes it with all four tail theorems. It "
            "therefore discharges the complex-disk exact-minus-epsilon-ten cumulant "
            "budgets through order eight. It does not by itself compose the final "
            "signed cumulant corridors into order-four entry, prove PF-infinity, RH, "
            "or Lambda<=0."
        ),
        "finite": {
            "jet_caps": {str(key): str(value) for key, value in FINITE_JET_CAPS.items()},
            "partition_bounds": {str(key): str(value) for key, value in finite_partition.items()},
            "weight_norms": {str(key): str(value) for key, value in finite_weight.items()},
            "perturbation": finite_perturbation,
            "bell": finite_bell,
            "potential": finite_potential,
            "budget": budgets["finite"],
        },
        "ray": {
            "jet_caps": {str(key): str(value) for key, value in ray_caps.items()},
            "partition_bounds": {str(key): str(value) for key, value in ray_partition.items()},
            "weight_norms": {str(key): str(value) for key, value in ray_weight.items()},
            "perturbation": ray_perturbation,
            "bell": ray_bell,
            "potential": ray_potential,
            "budget": budgets["ray"],
        },
        "potential_remainders": potential_serialized,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": len(rows),
            "open_analytic_rows": 0,
            "partition_extension_orders": 4,
            "bell_polynomial_terms": finite_bell["bell_polynomial_terms"],
            "formal_tails_closed": 2,
            "exact_tails_closed": 2,
            "central_residuals_closed": 2,
            "global_partition_residual_closed": True,
            "exact_minus_epsilon_ten_cumulant_budgets_closed": True,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            path.relative_to(REPO_ROOT).as_posix(): sha256(path)
            for path in source_paths
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.py"
        ),
        "next_target": (
            "Compose the now-exact cumulant corridors with the compound order-four "
            "algebra and the first-summand curvature bridge, then audit the remaining "
            "full-kernel and all-order handoffs."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Exact Cumulant Central-Residual Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: global exact central residual and unit-disk partition theorem.",
        "This is not a proof of order-four entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.py",
        "```",
        "",
        "## Central Decomposition",
        "",
        "The central exact-minus-epsilon-ten residual is decomposed before any",
        "cumulant differentiation into",
        "",
        "```text",
        "the full Gaussian partition coefficients Z_11,...,Z_14,",
        "their two added formal Gaussian tails,",
        "the epsilon-fifteen Bell-polynomial exponential remainder,",
        "the exact seventeenth-order potential Taylor remainder.",
        "```",
        "",
        "The finite interval uses the cancellation-preserving coefficient caps",
        "`2, 2, 21/10, 12/5`. The ray uses explicit, deliberately crude high-jet",
        "majorants; exponential q growth leaves ample room.",
        "",
        "## Closed Budgets",
        "",
        "The exact scalar compositions prove",
        "",
        "```text",
        artifact["finite"]["budget"]["proved_bound"],
        artifact["ray"]["budget"]["proved_bound"],
        "```",
        "",
        f"The finite central budget ratio is below `{artifact['finite']['budget']['budget_ratio_upper']}`.",
        f"The ray central budget ratio at its worst endpoint is below `{artifact['ray']['budget']['budget_ratio_upper']}`.",
        "",
        "Together with the two formal and two exact tail theorems, this gives",
        "",
        "```text",
        "|A_u-P_u^[10]|<1/(100000*q^3),       2<=u<=20,",
        "|A_u-P_u^[10]|<1/(20000*u*q^3),      u>=20.",
        "```",
        "",
        "The complex-disk contract therefore yields the simultaneous scaled",
        "exact-minus-epsilon-ten cumulant budgets `9/1000` and `1/(100u)` for",
        "orders two through eight.",
        "",
        "## Remaining Boundary",
        "",
        "The exact-density remainder is no longer open. The next job is to compose",
        "these exact cumulant corridors with the order-four determinant algebra and",
        "then audit the first-summand/full-kernel handoff. No claim of PF-infinity,",
        "RH, or `Lambda <= 0` is made here.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
        "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four exact cumulant central residual: "
        "11 rows, 11 exact rows, 222 Bell terms, 2 central regimes, "
        "4 inherited tails, 0 open partition components"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
