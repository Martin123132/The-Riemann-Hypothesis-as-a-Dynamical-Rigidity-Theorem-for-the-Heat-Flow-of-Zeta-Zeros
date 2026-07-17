#!/usr/bin/env python3
"""Prove the uniform signed order-six tail and reduce propagation to entry."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import itertools
import json
from math import factorial
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_uniform_tail_flow_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md"
)
ORDER5_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.json"
)
HEAT_TILT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.json"
)
HIGHER_THETA_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_superpolynomial_first_summand_dominance.json"
)
XI_RATIO_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.json"
)
ORDER = 6
BACKWARD_SHIFT_MAX = 2 * (ORDER - 1)
LEADING_DEGREE = ORDER * (ORDER - 1) // 2
MAX_RATIO_ORDER = LEADING_DEGREE + 1
RAW_LEADING_CONSTANT = -1132462080


@dataclass(frozen=True)
class ReductionRow:
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


def validate_sources() -> dict:
    order5 = load_json(ORDER5_SOURCE)
    heat_tilt = load_json(HEAT_TILT_SOURCE)
    higher_theta = load_json(HIGHER_THETA_SOURCE)
    xi_ratio = load_json(XI_RATIO_SOURCE)
    summary = order5.get("summary", {})
    if summary.get("contiguous_order_five_theorems") != 1:
        raise RuntimeError("strict order-five source is not closed")
    if summary.get("arbitrary_column_order_five_theorems") != 1:
        raise RuntimeError("arbitrary-column order-five source is not closed")
    if heat_tilt.get("summary", {}).get("uniform_heat_tilt_theorems") != 1:
        raise RuntimeError("published heat-tilt source is not closed")
    if higher_theta.get("summary", {}).get("open_analytic_rows") != 0:
        raise RuntimeError("uniform higher-theta source is not closed")
    if xi_ratio.get("eventual_theorem", {}).get("threshold_effective_here") is not False:
        raise RuntimeError("published Xi ratio source contract changed")
    return {
        "order5_status": order5.get("status"),
        "order5_sign": (
            "Q_(5,n)(lambda)=H_(5,n)(lambda)>0 for every n>=0 and "
            "-100<=lambda<=0"
        ),
        "heat_tilt_status": heat_tilt.get("status"),
        "higher_theta_status": higher_theta.get("status"),
        "xi_ratio_source": xi_ratio.get("published_ratio_input", {}).get("source"),
    }


def uniform_suitability_extension() -> dict:
    v, T = sp.symbols("v T")
    exponent = [sp.Integer(0)] * (MAX_RATIO_ORDER + 1)
    for degree in range(1, MAX_RATIO_ORDER + 1):
        log_coefficient = sp.Rational((-1) ** (degree + 1), degree)
        log_square_coefficient = (
            sp.Integer(0)
            if degree < 2
            else (-1) ** degree
            * sp.Rational(2, degree)
            * sp.harmonic(degree - 1)
        )
        exponent[degree] = sp.expand(
            -T * (v * log_coefficient / 8 + log_square_coefficient / 16)
        )
    coefficients = [sp.Integer(0)] * (MAX_RATIO_ORDER + 1)
    coefficients[0] = sp.Integer(1)
    rows = []
    for degree in range(1, MAX_RATIO_ORDER + 1):
        coefficients[degree] = sp.expand(
            sum(
                index * exponent[index] * coefficients[degree - index]
                for index in range(1, degree + 1)
            )
            / degree
        )
    for degree, coefficient in enumerate(coefficients):
        polynomial = sp.Poly(coefficient, T, v)
        if polynomial.degree(T) > degree or polynomial.degree(v) > degree:
            raise RuntimeError(f"heat suitability degree failed at {degree}")
        rows.append(
            {
                "degree": degree,
                "T_degree": polynomial.degree(T),
                "v_degree": polynomial.degree(v),
                "coefficient": str(sp.factor(coefficient)),
            }
        )
    return {
        "family": "f_T(t)=exp(-T*(log t)^2/16), 0<=T<=100",
        "uniform_bound": (
            "sup_(0<=T<=100)|f_(T,r)(v)|=O_r(v^r), r=0,...,16"
        ),
        "coefficient_checks": rows,
    }


def lambert_extension() -> dict:
    k, w = sp.symbols("k w", positive=True)

    def derivative(expression: sp.Expr) -> sp.Expr:
        return sp.factor(
            sp.diff(expression, k)
            + sp.diff(expression, w) * w / (k * (1 + w))
        )

    current = w**2
    rows = []
    for order in range(1, MAX_RATIO_ORDER + 1):
        current = derivative(current)
        normalized = sp.factor(current * k**order / w)
        limit = sp.limit(normalized, w, sp.oo)
        if not limit.is_finite or limit == 0:
            raise RuntimeError(f"Lambert derivative limit failed at {order}")
        rows.append(
            {
                "order": order,
                "derivative": str(current),
                "normalized": str(normalized),
                "limit_w_to_infinity": str(limit),
            }
        )
    return {
        "identity": "w'(k)=w/(k*(1+w)), w=W(2*k/pi)",
        "derivative_bound": (
            "d^m/dk^m w(k)^2=O(w(k)/k^m), m=1,...,16"
        ),
        "difference_bound": (
            "Delta_k^m log R_T^(1)(k)=O(log(k)/k^m), "
            "m=2,...,16, uniformly for 0<=T<=100"
        ),
        "rows": rows,
    }


def newton_transfer() -> dict:
    j = sp.symbols("j")
    differences = sp.symbols("D1:11")
    polynomial = sp.expand(
        sum(
            differences[order - 1] * sp.binomial(j, order).expand(func=True)
            for order in range(1, 11)
        )
    )
    rows = []
    for degree in range(1, 11):
        coefficient = sp.expand(polynomial).coeff(j, degree)
        support = [
            order
            for order, difference in enumerate(differences, start=1)
            if coefficient.coeff(difference) != 0
        ]
        diagonal = sp.factor(coefficient.coeff(differences[degree - 1]))
        if not support or min(support) != degree:
            raise RuntimeError(f"Newton support failed at degree {degree}")
        if diagonal != sp.Rational(1, sp.factorial(degree)):
            raise RuntimeError(f"Newton diagonal failed at degree {degree}")
        rows.append(
            {
                "degree": degree,
                "coefficient": str(sp.factor(coefficient)),
                "difference_orders": support,
                "diagonal": str(diagonal),
            }
        )
    return {
        "interpolant": "q_M(j)=sum_(m=1)^10 D_m(M)*binom(j,m), 0<=j<=10",
        "triangular_support": (
            "the coefficient of j^r uses only D_m with m>=r and diagonal D_r/r!"
        ),
        "heat_input": "D_m(M)=O(log(M)/M^m), m=2,...,10",
        "graded_output": (
            "[j^r]q_M(j)/h(M)^(r-1)=O(log(M)/M)=o(1), r=2,...,10"
        ),
        "coefficient_rows": rows,
    }


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def weighted_partitions(total: int, maximum_weight: int) -> list[tuple[int, ...]]:
    counts = [0] * maximum_weight
    rows: list[tuple[int, ...]] = []

    def visit(weight: int, remainder: int) -> None:
        if weight == 0:
            if remainder == 0:
                rows.append(tuple(counts))
            return
        for count in range(remainder // weight + 1):
            counts[weight - 1] = count
            visit(weight - 1, remainder - count * weight)
        counts[weight - 1] = 0

    visit(min(total, maximum_weight), total)
    return rows


def determinant_cancellation() -> dict:
    permutations = tuple(itertools.permutations(range(ORDER)))
    permutation_rows = []
    for permutation in permutations:
        power_sums = [
            sum(
                (BACKWARD_SHIFT_MAX - row - permutation[row]) ** (weight + 1)
                for row in range(ORDER)
            )
            for weight in range(1, LEADING_DEGREE + 1)
        ]
        permutation_rows.append((permutation_sign(permutation), power_sums))

    coefficient_layers: list[list[dict]] = []
    monomials_checked = 0
    for degree in range(LEADING_DEGREE + 1):
        nonzero = []
        for counts in weighted_partitions(degree, LEADING_DEGREE):
            monomials_checked += 1
            used = [
                (weight + 1, count)
                for weight, count in enumerate(counts)
                if count
            ]
            alternating_sum = 0
            for sign, power_sums in permutation_rows:
                term = sign
                for weight, count in used:
                    term *= power_sums[weight - 1] ** count
                alternating_sum += term
            if alternating_sum == 0:
                continue
            denominator = 1
            for count in counts:
                denominator *= factorial(count)
            coefficient = Fraction(
                (-1) ** sum(counts) * alternating_sum,
                denominator,
            )
            nonzero.append(
                {
                    "counts_G2_through_G16": list(counts),
                    "coefficient": str(coefficient),
                }
            )
        coefficient_layers.append(nonzero)

    expected_counts = [LEADING_DEGREE] + [0] * (LEADING_DEGREE - 1)
    expected = [
        {
            "counts_G2_through_G16": expected_counts,
            "coefficient": str(RAW_LEADING_CONSTANT),
        }
    ]
    if any(coefficient_layers[degree] for degree in range(LEADING_DEGREE)):
        raise RuntimeError("order-six determinant did not vanish below h^15")
    if coefficient_layers[LEADING_DEGREE] != expected:
        raise RuntimeError("unexpected order-six leading determinant coefficient")

    vandermonde_constant = (
        (-1) ** LEADING_DEGREE
        * 2**LEADING_DEGREE
        * sp.prod(sp.factorial(index) for index in range(1, ORDER))
    )
    if vandermonde_constant != RAW_LEADING_CONSTANT:
        raise RuntimeError("order-six Vandermonde constant changed")
    return {
        "matrix_entry": (
            "K_(i,j)(h)=exp(-sum_(m>=2)G_m*h^(m-1)*(10-i-j)^m)"
        ),
        "truncation": "m=2,...,16 and powers h^0,...,h^15",
        "coefficients_h0_through_h15": ["0"] * LEADING_DEGREE
        + ["-1132462080*G2**15"],
        "first_nonzero_term": "-1132462080*G2^15*h^15",
        "signed_first_nonzero_term": "1132462080*G2^15*h^15 for Q_6=-H_6",
        "general_constant": (
            "epsilon_6*2^binom(6,2)*prod_(j=1)^5 j!=-1132462080"
        ),
        "independence": "G3,...,G16 cancel from the first nonzero coefficient",
        "permutations_checked": len(permutations),
        "weighted_monomials_checked": monomials_checked,
        "permutation_monomial_checks": len(permutations) * monomials_checked,
        "nonzero_layers": coefficient_layers,
    }


def determinant_polynomial(
    values: tuple[sp.Symbol, ...], order: int, shift: int
) -> sp.Expr:
    return sp.Matrix(
        [[values[shift + i + j] for j in range(order)] for i in range(order)]
    ).det(method="domain-ge")


def delta_derivative(expression: sp.Expr, values: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(expression, values[index]) * values[index + 1]
            for index in range(len(values) - 1)
        )
    )


def exact_flow_and_condensation() -> dict:
    values = sp.symbols("a0:13")
    h4_n2 = determinant_polynomial(values, 4, 2)
    h5_n0 = determinant_polynomial(values, 5, 0)
    h5_n1 = determinant_polynomial(values, 5, 1)
    h5_n2 = determinant_polynomial(values, 5, 2)
    h6_n0 = determinant_polynomial(values, 6, 0)
    h6_n1 = determinant_polynomial(values, 6, 1)
    condensation = sp.cancel(
        h6_n0 * h4_n2 - (h5_n0 * h5_n2 - h5_n1**2)
    )
    if condensation != 0:
        raise RuntimeError("order-six condensation identity failed")
    plucker = sp.cancel(
        h5_n1 * delta_derivative(h6_n0, values)
        - h5_n0 * h6_n1
        - delta_derivative(h5_n1, values) * h6_n0
    )
    if plucker != 0:
        raise RuntimeError("order-six Plucker flow identity failed")
    heat_derivative = sp.expand(
        sum(
            sp.diff(h6_n0, values[index])
            * (4 * index + 2)
            * values[index + 1]
            for index in range(len(values) - 1)
        )
    )
    affine = sp.cancel(heat_derivative - 42 * delta_derivative(h6_n0, values))
    if affine != 0:
        raise RuntimeError("order-six affine heat derivative failed")
    return {
        "signed_coordinate": "Q_(m,n)=epsilon_m*H_(m,n), epsilon_m=(-1)^binom(m,2)",
        "signed_condensation_recursion": (
            "Q_(m,n)*Q_(m-2,n+2)=Q_(m-1,n+1)^2-"
            "Q_(m-1,n)*Q_(m-1,n+2)"
        ),
        "general_entry_equivalence": (
            "Q_(m,n)>0 iff Q_(m-1,n+1)^2>"
            "Q_(m-1,n)*Q_(m-1,n+2), provided Q_(m-2,n+2)>0"
        ),
        "general_affine_derivative": (
            "Q_(m,n)'=c_(m,n)*delta(Q_(m,n)), c_(m,n)=4*n+8*m-6"
        ),
        "general_plucker": (
            "Q_(m-1,n+1)*delta(Q_(m,n))=Q_(m-1,n)*Q_(m,n+1)+"
            "delta(Q_(m-1,n+1))*Q_(m,n)"
        ),
        "general_cooperative_flow": (
            "Q_(m,n)'=a_(m,n)*Q_(m,n+1)+b_(m,n)*Q_(m,n), "
            "a_(m,n)=c_(m,n)*Q_(m-1,n)/Q_(m-1,n+1)>0, "
            "b_(m,n)=c_(m,n)/(c_(m,n)-4)*(log Q_(m-1,n+1))'"
        ),
        "order6_condensation": (
            "H_(6,n)*H_(4,n+2)=H_(5,n)*H_(5,n+2)-H_(5,n+1)^2"
        ),
        "order6_entry_equivalence": (
            "Q_(6,n)=-H_(6,n)>0 iff H_(5,n+1)^2>"
            "H_(5,n)*H_(5,n+2)"
        ),
        "order6_affine_derivative": "Q_(6,n)'=(4*n+42)*delta(Q_(6,n))",
        "order6_cooperative_flow": (
            "Q_n'=a_n*Q_(n+1)+b_n*Q_n, "
            "a_n=(4*n+42)*H_(5,n)/H_(5,n+1)>0, "
            "b_n=((4*n+42)/(4*n+38))*(log H_(5,n+1))'"
        ),
        "symbolic_residuals": {
            "condensation": str(condensation),
            "plucker": str(plucker),
            "affine_heat": str(affine),
        },
    }


def lower_layer_countermodel() -> dict:
    values = tuple(
        [sp.Rational(1, sp.factorial(index)) for index in range(10)]
        + [sp.Rational(1, 3590000)]
    )

    def determinant(order: int, shift: int) -> sp.Rational:
        return sp.factor(
            sp.Matrix(
                [
                    [values[shift + i + j] for j in range(order)]
                    for i in range(order)
                ]
            ).det()
        )

    signed_rows = {}
    for order in range(1, 6):
        orientation = (-1) ** (order * (order - 1) // 2)
        count = len(values) - 2 * order + 2
        rows = [
            sp.factor(orientation * determinant(order, shift))
            for shift in range(count)
        ]
        if not all(value > 0 for value in rows):
            raise RuntimeError(f"countermodel lower layer failed at order {order}")
        signed_rows[str(order)] = [str(value) for value in rows]
    h6 = determinant(6, 0)
    q6 = -h6
    h5_margin = sp.factor(
        determinant(5, 1) ** 2 - determinant(5, 0) * determinant(5, 2)
    )
    if h6 <= 0 or q6 >= 0 or h5_margin >= 0:
        raise RuntimeError("countermodel did not violate signed order six")
    return {
        "sequence": [str(value) for value in values],
        "strict_signed_lower_layers": signed_rows,
        "H6_n0": str(h6),
        "Q6_n0": str(q6),
        "H5_log_concavity_margin": str(h5_margin),
        "conclusion": (
            "strict signed contiguous layers through order five do not imply order six"
        ),
    }


def conclusions(cancellation: dict, flow: dict) -> dict:
    if cancellation.get("first_nonzero_term") != (
        "-1132462080*G2^15*h^15"
    ):
        raise RuntimeError("order-six main term changed")
    return {
        "heat_tilt_extension": (
            "Delta_k^m log(A_k^(1)(-T)/A_k^(1)(0))=O(log(k)/k^m), "
            "m=2,...,16, uniformly for 0<=T<=100"
        ),
        "higher_theta": (
            "all fixed local differences of log(A_k(-T)/A_k^(1)(-T)) "
            "are O_p,m(k^-p) uniformly"
        ),
        "uniform_ratio": (
            "the complete family has the compact-uniform order-six graded ratio "
            "contract with G_2(T,M)->1"
        ),
        "determinant_asymptotic": (
            "Q_(6,n)(lambda)=positive_scale(lambda,n)*(1132462080*"
            "G_2(lambda,n)^15*h(lambda,n)^15+o(h(lambda,n)^15)) uniformly"
        ),
        "uniform_eventual_tail": (
            "there exists N_6 such that Q_(6,n)(lambda)=-H_(6,n)(lambda)>0 "
            "for every n>=N_6 and -100<=lambda<=0"
        ),
        "conditional_entry": (
            "Q_(6,n)(-100)>0 for every integer n>=0, equivalently "
            "H_(5,n+1)(-100)^2>H_(5,n)(-100)*H_(5,n+2)(-100)"
        ),
        "variation_of_constants": (
            "Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda "
            "E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)"
        ),
        "conditional_forward": (
            "[Q_(6,n)(-100)>0 for all n] => [Q_(6,n)(lambda)>0 for all n "
            "and -100<=lambda<=0]"
        ),
        "conditional_arbitrary_columns": (
            "completed signed contiguous layers through order six imply every "
            "arbitrary-column signed layer through order six"
        ),
        "open_entry_target": flow["order6_entry_equivalence"] + (
            " at lambda=-100 for every n>=0"
        ),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    suitability = uniform_suitability_extension()
    lambert = lambert_extension()
    newton = newton_transfer()
    cancellation = determinant_cancellation()
    flow = exact_flow_and_condensation()
    countermodel = lower_layer_countermodel()
    exact = conclusions(cancellation, flow)
    rows = [
        ReductionRow(
            "co6utfr_01_order5_input",
            "theorem_input",
            "ready_to_apply",
            "The complete signed layer through fixed order five supplies the positive lower cone for order six.",
            sources["order5_sign"],
            "Completed fixed-order input only; no order-six sign is imported.",
        ),
        ReductionRow(
            "co6utfr_02_heat_tilt_extension",
            "published_theorem_composition",
            "ready_to_apply",
            "The all-order suitable-multiplier theorem extends the compact heat-tilt audit through order sixteen.",
            exact["heat_tilt_extension"],
            "First theta summand only; compact heat interval.",
            {"suitability": suitability, "lambert": lambert},
        ),
        ReductionRow(
            "co6utfr_03_higher_theta",
            "theorem_input",
            "ready_to_apply",
            "Uniform superpolynomial first-summand dominance removes every higher-theta correction needed by the order-six stencil.",
            exact["higher_theta"],
            "Complete-to-first-summand correction only.",
        ),
        ReductionRow(
            "co6utfr_04_newton_transfer",
            "exact_symbolic_lemma",
            "ready_to_apply",
            "Eleven-node Newton interpolation transfers local heat differences to the graded order-six ratio coefficients.",
            newton["graded_output"],
            "Exact on backward shifts zero through ten.",
            newton,
        ),
        ReductionRow(
            "co6utfr_05_determinant_cancellation",
            "exact_symbolic_theorem",
            "ready_to_apply",
            "The normalized six-by-six determinant vanishes through h^14 and has the universal signed first coefficient.",
            cancellation["first_nonzero_term"],
            "Exact 720-permutation, 684-weighted-monomial audit.",
            cancellation,
        ),
        ReductionRow(
            "co6utfr_06_uniform_eventual_tail",
            "exact_theorem_composition",
            "ready_to_apply",
            "The heat-deformed Xi family has one eventual positive signed order-six tail uniformly through lambda zero.",
            exact["uniform_eventual_tail"],
            "Non-effective eventual signed theorem only.",
            {"asymptotic": exact["determinant_asymptotic"]},
        ),
        ReductionRow(
            "co6utfr_07_signed_condensation_recursion",
            "exact_identity",
            "ready_to_apply",
            "Every new signed contiguous layer is strict log-concavity of the preceding signed layer.",
            flow["signed_condensation_recursion"] + "; " + flow["general_entry_equivalence"],
            "General Desnanot-Jacobi identity with exact orientation parity.",
        ),
        ReductionRow(
            "co6utfr_08_order6_entry_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Signed order-six entry is exactly strict log-concavity of the completed positive H_5 sequence.",
            flow["order6_condensation"] + "; " + flow["order6_entry_equivalence"],
            "Exact endpoint coordinate; no positivity asserted here.",
        ),
        ReductionRow(
            "co6utfr_09_general_cooperative_recursion",
            "exact_identity",
            "ready_to_apply",
            "At every fixed order, the heat flow is cooperative over the completed preceding signed layer and needs no next-order sign.",
            flow["general_cooperative_flow"],
            "General affine-Hankel and adjacent Plucker algebra.",
        ),
        ReductionRow(
            "co6utfr_10_order6_cooperative_flow",
            "exact_forward_lemma",
            "ready_to_apply",
            "The signed order-six heat flow has a strictly positive nearest-neighbor coefficient throughout the completed order-five cone.",
            flow["order6_cooperative_flow"],
            "Exact order-six specialization with three zero symbolic residuals.",
            flow["symbolic_residuals"],
        ),
        ReductionRow(
            "co6utfr_11_conditional_forward",
            "exact_conditional_theorem",
            "conditional_on_open_input",
            "All-shift lambda=-100 entry would complete signed contiguous and arbitrary-column order six through lambda zero.",
            exact["conditional_forward"] + "; " + exact["conditional_arbitrary_columns"],
            "Conditional on the explicit endpoint theorem only.",
        ),
        ReductionRow(
            "co6utfr_12_countermodel",
            "countermodel_gate",
            "ready_to_apply",
            "The completed lower signed layers do not imply order six abstractly.",
            "Q_(6,0)=" + countermodel["Q6_n0"] + "<0",
            "Exact rational finite sequence; blocks lower-layer promotion.",
            countermodel,
        ),
        ReductionRow(
            "co6utfr_13_open_entry",
            "analytic_theorem_target",
            "not_ready_to_apply",
            "Prove strict all-shift signed order-six entry at lambda=-100.",
            exact["conditional_entry"],
            "Sole open input in this fixed-order-six propagation route.",
        ),
    ]
    source_paths = (ORDER5_SOURCE, HEAT_TILT_SOURCE, HIGHER_THETA_SOURCE, XI_RATIO_SOURCE)
    return {
        "kind": "jensen_window_pf_compound_order6_uniform_tail_flow_reduction",
        "date": "2026-07-13",
        "status": (
            "exact uniform eventual signed order-six tail and conditional forward "
            "reduction with one open lambda=-100 entry"
        ),
        "proof_boundary": (
            "This artifact proves the signed order-six uniform eventual tail, the "
            "general condensation/cooperative recursion, and a conditional forward "
            "theorem. It does not prove all-shift order-six entry, full order-six "
            "invariance, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [path.relative_to(REPO_ROOT).as_posix() for path in source_paths],
        "source_contract": sources,
        "uniform_suitability_extension": suitability,
        "lambert_extension": lambert,
        "newton_transfer": newton,
        "determinant_cancellation": cancellation,
        "exact_flow": flow,
        "countermodel": countermodel,
        "conclusions": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": 11,
            "conditional_rows": 1,
            "open_rows": 1,
            "suitability_coefficients_checked": len(suitability["coefficient_checks"]),
            "lambert_derivative_orders_checked": len(lambert["rows"]),
            "newton_coefficients_checked": len(newton["coefficient_rows"]),
            "determinant_permutations_checked": cancellation["permutations_checked"],
            "weighted_monomials_checked": cancellation["weighted_monomials_checked"],
            "uniform_eventual_tail_theorems": 1,
            "signed_condensation_recursions": 1,
            "cooperative_flow_recursions": 1,
            "conditional_forward_theorems": 1,
            "lower_layer_countermodels": 1,
            "open_entry_targets": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order6_uniform_tail_flow_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order6_uniform_tail_flow_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    cancellation = artifact["determinant_cancellation"]
    flow = artifact["exact_flow"]
    exact = artifact["conclusions"]
    countermodel = artifact["countermodel"]
    lines = [
        "# Jensen-Window PF Compound Order-Six Uniform Tail And Flow Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact uniform eventual signed order-six tail and conditional",
        "forward reduction with one open lambda=-100 entry. This is not a proof",
        "of all-shift order six, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_uniform_tail_flow_reduction.py",
        "```",
        "",
        "## Universal Signed Tail",
        "",
        "The published all-order Xi ratio expansion, a compact heat-tilt audit",
        "through order sixteen, and superpolynomial higher-theta suppression",
        "reduce the normalized six-by-six determinant to",
        "",
        "```text",
        cancellation["matrix_entry"],
        "[h^0,...,h^15] det K=[0,...,0,-1132462080*G_2^15].",
        "```",
        "",
        f"The exact audit checks `{cancellation['permutations_checked']}` permutations",
        f"against `{cancellation['weighted_monomials_checked']}` weighted monomials,",
        f"or `{cancellation['permutation_monomial_checks']}` exact permutation-monomial",
        "terms. Every contribution containing `G_3,...,G_16` cancels from the",
        "first nonzero coefficient. Since `Q_6=-H_6`, this proves",
        "",
        "```text",
        exact["determinant_asymptotic"],
        exact["uniform_eventual_tail"],
        "```",
        "",
        "The threshold is finite but non-effective.",
        "",
        "## Fixed-Order Recursion",
        "",
        "Put `epsilon_m=(-1)^binom(m,2)` and `Q_(m,n)=epsilon_m*H_(m,n)`.",
        "Desnanot-Jacobi and the orientation identity",
        "`epsilon_m*epsilon_(m-2)=-1` give",
        "",
        "```text",
        flow["signed_condensation_recursion"],
        flow["general_entry_equivalence"],
        "```",
        "",
        "Thus every new signed contiguous layer is strict log-concavity of the",
        "preceding signed layer. The affine heat derivative and adjacent Plucker",
        "identity similarly give the general cooperative recursion",
        "",
        "```text",
        flow["general_affine_derivative"],
        flow["general_plucker"],
        flow["general_cooperative_flow"],
        "```",
        "",
        "The positive off-diagonal coefficient uses only the completed order",
        "`m-1` sign; no order `m+1` sign is required.",
        "",
        "## Order-Six Specialization",
        "",
        "At order six the endpoint coordinate and flow are",
        "",
        "```text",
        flow["order6_condensation"],
        flow["order6_entry_equivalence"],
        flow["order6_cooperative_flow"],
        "```",
        "",
        "The completed `H_5>0` layer makes `a_n>0`. The uniform signed tail",
        "then confines a hypothetical first loss to finitely many shifts, so",
        "variation of constants proves the conditional theorem",
        "",
        "```text",
        exact["conditional_forward"],
        exact["conditional_arbitrary_columns"],
        "```",
        "",
        "## Countermodel Gate",
        "",
        "The rational sequence",
        "",
        "```text",
        "(" + ",".join(countermodel["sequence"]) + ")",
        "```",
        "",
        "has every available strict signed contiguous minor through order five,",
        "but",
        "",
        "```text",
        "H_(6,0)=" + countermodel["H6_n0"] + ">0,",
        "Q_(6,0)=" + countermodel["Q6_n0"] + "<0,",
        "H_(5,1)^2-H_(5,0)H_(5,2)=" + countermodel["H5_log_concavity_margin"] + "<0.",
        "```",
        "",
        "Order six is therefore genuinely new and cannot be promoted from the",
        "completed lower cone.",
        "",
        "## Remaining Endpoint Target",
        "",
        "The remaining fixed-order problem is all-shift signed order-six entry:",
        "",
        "```text",
        exact["conditional_entry"],
        "```",
        "",
        "A separate rigorous prefix certificate handles the currently available",
        "finite coefficient range. The analytic all-shift tail remains open.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md",
        "outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md",
        "outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-six uniform tail and flow reduction: "
        f"{summary['rows']} rows, "
        f"{summary['determinant_permutations_checked']} permutations, "
        f"{summary['weighted_monomials_checked']} weighted monomials, "
        f"{summary['open_entry_targets']} open entry"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
