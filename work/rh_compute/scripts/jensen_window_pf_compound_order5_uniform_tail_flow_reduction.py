#!/usr/bin/env python3
"""Prove the uniform order-five tail and reduce propagation to lambda=-100 entry."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import itertools
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_uniform_tail_flow_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md"
)
ORDER4_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.json"
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
ASYMPTOTIC_ORDER = 10


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
    order4 = load_json(ORDER4_SOURCE)
    heat_tilt = load_json(HEAT_TILT_SOURCE)
    higher_theta = load_json(HIGHER_THETA_SOURCE)
    xi_ratio = load_json(XI_RATIO_SOURCE)
    if order4.get("summary", {}).get("lambda_zero_all_shift_theorems") != 1:
        raise RuntimeError("strict order-four source is not closed")
    if heat_tilt.get("summary", {}).get("uniform_heat_tilt_theorems") != 1:
        raise RuntimeError("uniform heat-tilt source is not closed")
    if higher_theta.get("summary", {}).get("open_analytic_rows") != 0:
        raise RuntimeError("uniform higher-theta source is not closed")
    if xi_ratio.get("eventual_theorem", {}).get("threshold_effective_here") is not False:
        raise RuntimeError("published Xi ratio source contract changed")
    return {
        "order4_status": order4.get("status"),
        "heat_tilt_status": heat_tilt.get("status"),
        "higher_theta_status": higher_theta.get("status"),
        "xi_ratio_source": xi_ratio.get("published_ratio_input", {}).get("source"),
        "order4_sign": "H_(4,n)(lambda)>0 for every n>=0 and -100<=lambda<=0",
    }


def uniform_suitability_extension() -> dict:
    v, T = sp.symbols("v T")
    exponent_coefficients = [sp.Integer(0)] * 12
    for degree in range(1, 12):
        log_coefficient = sp.Rational((-1) ** (degree + 1), degree)
        log_square_coefficient = (
            sp.Integer(0)
            if degree < 2
            else (-1) ** degree
            * sp.Rational(2, degree)
            * sp.harmonic(degree - 1)
        )
        exponent_coefficients[degree] = sp.expand(
            -T * (v * log_coefficient / 8 + log_square_coefficient / 16)
        )
    coefficients = [sp.Integer(0)] * 12
    coefficients[0] = sp.Integer(1)
    for degree in range(1, 12):
        coefficients[degree] = sp.expand(
            sum(
                k * exponent_coefficients[k] * coefficients[degree - k]
                for k in range(1, degree + 1)
            )
            / degree
        )
    checks = []
    for degree in range(12):
        coefficient = sp.factor(coefficients[degree])
        polynomial = sp.Poly(coefficient, T, v)
        if polynomial.degree(T) > degree or polynomial.degree(v) > degree:
            raise RuntimeError(f"heat suitability degree failed at {degree}")
        checks.append(
            {
                "degree": degree,
                "T_degree": polynomial.degree(T),
                "v_degree": polynomial.degree(v),
                "coefficient": str(coefficient),
            }
        )
    return {
        "family": "f_T(t)=exp(-T*(log t)^2/16), 0<=T<=100",
        "ratio_identity": (
            "f_T(t*(1+x))/f_T(t)=exp(-T*((v/8)*log(1+x)+"
            "log(1+x)^2/16)), t=exp(v)"
        ),
        "uniform_bound": (
            "sup_(0<=T<=100)|f_(T,r)(v)|=O_r(v^r), r=0,...,11"
        ),
        "coefficient_checks": checks,
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
    for order in range(1, 12):
        current = derivative(current)
        normalized = sp.factor(current * k**order / w)
        limit = sp.limit(normalized, w, sp.oo)
        if not limit.is_finite:
            raise RuntimeError(f"Lambert derivative limit failed at order {order}")
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
            "d^m/dk^m w(k)^2=O(w(k)/k^m), m=1,...,11"
        ),
        "difference_bound": (
            "Delta_k^m log R_T^(1)(k)=O(log(k)/k^m), "
            "m=2,...,11, uniformly for 0<=T<=100"
        ),
        "rows": rows,
    }


def newton_transfer() -> dict:
    j = sp.symbols("j")
    differences = sp.symbols("D1:9")
    polynomial = sp.expand(
        sum(
            differences[order - 1] * sp.binomial(j, order).expand(func=True)
            for order in range(1, 9)
        )
    )
    rows = []
    for degree in range(1, 9):
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
        "interpolant": "q_M(j)=sum_(m=1)^8 D_m(M)*binom(j,m), 0<=j<=8",
        "triangular_support": (
            "the coefficient of j^r uses only D_m with m>=r and diagonal D_r/r!"
        ),
        "heat_input": "D_m(M)=O(log(M)/M^m), m=2,...,8",
        "graded_output": (
            "[j^r]q_M(j)/h(M)^(r-1)=O(log(M)/M)=o(1), r=2,...,8"
        ),
        "polynomial": str(polynomial),
        "coefficient_rows": rows,
    }


def series_multiply(left: list[sp.Expr], right: list[sp.Expr]) -> list[sp.Expr]:
    return [
        sp.expand(sum(left[k] * right[n - k] for k in range(n + 1)))
        for n in range(ASYMPTOTIC_ORDER + 1)
    ]


def exp_log_series(backward_shift: int, symbols: tuple[sp.Symbol, ...]) -> list[sp.Expr]:
    log_coefficients = [sp.Integer(0)] * (ASYMPTOTIC_ORDER + 1)
    for order in range(2, ASYMPTOTIC_ORDER + 2):
        log_coefficients[order - 1] = -(
            symbols[order - 2] * sp.Integer(backward_shift) ** order
        )
    exponential = [sp.Integer(0)] * (ASYMPTOTIC_ORDER + 1)
    exponential[0] = sp.Integer(1)
    for degree in range(1, ASYMPTOTIC_ORDER + 1):
        exponential[degree] = sp.expand(
            sum(
                k * log_coefficients[k] * exponential[degree - k]
                for k in range(1, degree + 1)
            )
            / degree
        )
    return exponential


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant_cancellation() -> dict:
    symbols = sp.symbols("G2:12")
    entry_series = {
        backward_shift: exp_log_series(backward_shift, symbols)
        for backward_shift in range(9)
    }
    determinant = [sp.Integer(0)] * (ASYMPTOTIC_ORDER + 1)
    permutations = tuple(itertools.permutations(range(5)))
    for permutation in permutations:
        product = [sp.Integer(1)] + [sp.Integer(0)] * ASYMPTOTIC_ORDER
        for row in range(5):
            backward_shift = 8 - row - permutation[row]
            product = series_multiply(product, entry_series[backward_shift])
        sign = permutation_sign(permutation)
        determinant = [
            sp.expand(current + sign * contribution)
            for current, contribution in zip(determinant, product)
        ]
    factored = [sp.factor(value) for value in determinant]
    expected = [sp.Integer(0)] * ASYMPTOTIC_ORDER + [
        294912 * symbols[0] ** 10
    ]
    if factored != expected:
        raise RuntimeError(f"unexpected order-five determinant expansion: {factored}")
    vandermonde_constant = (
        2 ** 10 * sp.prod(sp.factorial(j) for j in range(1, 5))
    )
    if vandermonde_constant != 294912:
        raise RuntimeError("order-five Vandermonde constant changed")
    return {
        "matrix_entry": (
            "K_(i,j)(h)=exp(-sum_(m>=2)G_m*h^(m-1)*(8-i-j)^m)"
        ),
        "truncation": "m=2,...,11 and powers h^0,...,h^10",
        "coefficients_h0_through_h10": [str(value) for value in factored],
        "first_nonzero_term": "294912*G2^10*h^10",
        "general_constant": (
            "epsilon_5*2^binom(5,2)*prod_(j=1)^4 j!=294912"
        ),
        "independence": "G3,...,G11 cancel from the first nonzero coefficient",
        "permutations_checked": len(permutations),
    }


def determinant_polynomial(values: tuple[sp.Symbol, ...], order: int, shift: int) -> sp.Expr:
    return sp.det(
        sp.Matrix(
            [[values[shift + i + j] for j in range(order)] for i in range(order)]
        )
    )


def delta_derivative(expression: sp.Expr, values: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(expression, values[index]) * values[index + 1]
            for index in range(len(values) - 1)
        )
    )


def exact_flow_and_condensation() -> dict:
    values = sp.symbols("a0:12")
    h3_n2 = determinant_polynomial(values, 3, 2)
    h4_n0 = determinant_polynomial(values, 4, 0)
    h4_n1 = determinant_polynomial(values, 4, 1)
    h4_n2 = determinant_polynomial(values, 4, 2)
    h5_n0 = determinant_polynomial(values, 5, 0)
    h5_n1 = determinant_polynomial(values, 5, 1)
    condensation = sp.factor(
        h5_n0 * h3_n2 - (h4_n0 * h4_n2 - h4_n1**2)
    )
    if condensation != 0:
        raise RuntimeError("order-five condensation identity failed")
    plucker = sp.factor(
        h4_n1 * delta_derivative(h5_n0, values)
        - h4_n0 * h5_n1
        - delta_derivative(h4_n1, values) * h5_n0
    )
    if plucker != 0:
        raise RuntimeError("order-five Plucker flow identity failed")
    heat_derivative = sp.expand(
        sum(
            sp.diff(h5_n0, values[index])
            * (4 * index + 2)
            * values[index + 1]
            for index in range(len(values) - 1)
        )
    )
    affine_residual = sp.factor(
        heat_derivative - 34 * delta_derivative(h5_n0, values)
    )
    if affine_residual != 0:
        raise RuntimeError("order-five affine heat derivative failed")
    return {
        "condensation": (
            "H_(5,n)*H_(3,n+2)=H_(4,n)*H_(4,n+2)-H_(4,n+1)^2"
        ),
        "entry_equivalence": (
            "because H_(3,n+2)<0 and H_(4,j)>0, H_(5,n)>0 iff "
            "H_(4,n+1)^2>H_(4,n)*H_(4,n+2)"
        ),
        "affine_derivative": "H_(5,n)'=(4*n+34)*delta(H_(5,n))",
        "plucker": (
            "H_(4,n+1)*delta(H_(5,n))=H_(4,n)*H_(5,n+1)+"
            "delta(H_(4,n+1))*H_(5,n)"
        ),
        "cooperative_flow": (
            "Q_n'=a_n*Q_(n+1)+b_n*Q_n, "
            "a_n=(4*n+34)*H_(4,n)/H_(4,n+1)>0, "
            "b_n=((4*n+34)/(4*n+30))*(log H_(4,n+1))'"
        ),
        "symbolic_residuals": {
            "condensation": str(condensation),
            "plucker": str(plucker),
            "affine_heat": str(affine_residual),
        },
    }


def lower_layer_countermodel() -> dict:
    values = tuple(
        [sp.Rational(1, sp.factorial(k)) for k in range(8)]
        + [sp.Rational(1, 42000)]
    )

    def determinant(order: int, shift: int) -> sp.Rational:
        return sp.det(
            sp.Matrix(
                [
                    [values[shift + i + j] for j in range(order)]
                    for i in range(order)
                ]
            )
        )

    signed_rows = {}
    for order in range(1, 5):
        orientation = (-1) ** (order * (order - 1) // 2)
        count = len(values) - (2 * order - 1) + 1
        rows = [sp.factor(orientation * determinant(order, n)) for n in range(count)]
        if not all(value > 0 for value in rows):
            raise RuntimeError(f"countermodel lower layer failed at order {order}")
        signed_rows[str(order)] = [str(value) for value in rows]
    h5 = sp.factor(determinant(5, 0))
    margin = sp.factor(determinant(4, 1) ** 2 - determinant(4, 0) * determinant(4, 2))
    if h5 >= 0 or margin >= 0:
        raise RuntimeError("countermodel did not violate order five")
    return {
        "sequence": [str(value) for value in values],
        "strict_signed_lower_layers": signed_rows,
        "H5_n0": str(h5),
        "H4_log_concavity_margin": str(margin),
        "conclusion": (
            "strict signed contiguous layers through order four do not imply order five"
        ),
    }


def asymptotic_and_flow_conclusions(cancellation: dict, flow: dict) -> dict:
    if cancellation["first_nonzero_term"] != "294912*G2^10*h^10":
        raise RuntimeError("order-five main term changed")
    return {
        "heat_tilt_extension": (
            "Delta_k^m log(A_k^(1)(-T)/A_k^(1)(0))=O(log(k)/k^m), "
            "m=2,...,11, uniformly for 0<=T<=100"
        ),
        "higher_theta": (
            "all fixed local differences of log(A_k(-T)/A_k^(1)(-T)) "
            "are O_p,m(k^-p) uniformly"
        ),
        "uniform_ratio": (
            "the complete family has the compact-uniform order-five graded "
            "ratio contract with G_2(T,M)->1"
        ),
        "determinant_asymptotic": (
            "H_(5,n)(lambda)=positive_scale(lambda,n)*(294912*"
            "G_2(lambda,n)^10*h(lambda,n)^10+o(h(lambda,n)^10)) uniformly"
        ),
        "uniform_eventual_tail": (
            "there exists N_5 such that H_(5,n)(lambda)>0 for every n>=N_5 "
            "and -100<=lambda<=0"
        ),
        "conditional_entry": "H_(5,n)(-100)>0 for every integer n>=0",
        "variation_of_constants": (
            "Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda "
            "E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)"
        ),
        "conditional_forward": (
            "[H_(5,n)(-100)>0 for all n] => [H_(5,n)(lambda)>0 for all n "
            "and -100<=lambda<=0]"
        ),
        "open_entry_target": flow["entry_equivalence"] + " at lambda=-100 for every n>=0",
    }


def build_artifact() -> dict:
    sources = validate_sources()
    suitability = uniform_suitability_extension()
    lambert = lambert_extension()
    newton = newton_transfer()
    cancellation = determinant_cancellation()
    flow = exact_flow_and_condensation()
    countermodel = lower_layer_countermodel()
    conclusions = asymptotic_and_flow_conclusions(cancellation, flow)
    rows = [
        ReductionRow(
            id="co5utfr_01_published_ratio_input",
            role="published_theorem_input",
            readiness="ready_to_apply",
            claim="The published Xi coefficient-ratio theorem supplies the all-order lambda-zero graded expansion.",
            formula="log(gamma(M-j)/gamma(M))=-sum_(m>=1)G_m(M)*Delta(M)^(2m-2)*j^m",
            proof_boundary="Published fixed-shift Xi asymptotics; no RH input.",
        ),
        ReductionRow(
            id="co5utfr_02_heat_tilt_extension",
            role="published_theorem_composition",
            readiness="ready_to_apply",
            claim="O'Sullivan's all-order suitable-multiplier theorem extends the compact heat-tilt differences through order eleven.",
            formula=conclusions["heat_tilt_extension"],
            proof_boundary="First theta summand only; compact heat interval.",
            diagnostics={"suitability": suitability, "lambert": lambert},
        ),
        ReductionRow(
            id="co5utfr_03_higher_theta",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="Uniform superpolynomial first-summand dominance removes all higher theta summands to every order needed here.",
            formula=conclusions["higher_theta"],
            proof_boundary="Complete-to-first-summand ratio correction only.",
        ),
        ReductionRow(
            id="co5utfr_04_newton_transfer",
            role="exact_symbolic_lemma",
            readiness="ready_to_apply",
            claim="Nine-node Newton interpolation transfers the heat differences to the graded order-five ratio coefficients.",
            formula=newton["graded_output"],
            proof_boundary="Exact on the nine determinant nodes.",
            diagnostics=newton,
        ),
        ReductionRow(
            id="co5utfr_05_determinant_cancellation",
            role="exact_symbolic_theorem",
            readiness="ready_to_apply",
            claim="The normalized five-by-five determinant vanishes through h^9 and has a universal positive first coefficient.",
            formula=cancellation["first_nonzero_term"],
            proof_boundary="Exact 120-permutation expansion through the first nonzero order.",
            diagnostics=cancellation,
        ),
        ReductionRow(
            id="co5utfr_06_uniform_eventual_tail",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The actual heat-deformed Xi family has one eventual positive contiguous order-five tail uniformly through lambda zero.",
            formula=conclusions["uniform_eventual_tail"],
            proof_boundary="Non-effective eventual theorem only.",
            diagnostics={"asymptotic": conclusions["determinant_asymptotic"]},
        ),
        ReductionRow(
            id="co5utfr_07_condensation",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Order-five entry is exactly strict log-concavity of the completed positive order-four determinants.",
            formula=flow["condensation"] + "; " + flow["entry_equivalence"],
            proof_boundary="Exact contiguous Desnanot-Jacobi identity.",
        ),
        ReductionRow(
            id="co5utfr_08_cooperative_flow",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The order-five heat flow is cooperative over the completed positive order-four layer and needs no order-six sign.",
            formula=flow["cooperative_flow"],
            proof_boundary="Exact affine-Hankel and Plucker algebra.",
            diagnostics=flow,
        ),
        ReductionRow(
            id="co5utfr_09_finite_confinement",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The uniform positive tail reduces forward propagation to finite backward induction from its boundary.",
            formula=conclusions["variation_of_constants"],
            proof_boundary="Requires strict lambda=-100 entry at the finitely many remaining shifts.",
        ),
        ReductionRow(
            id="co5utfr_10_conditional_forward",
            role="conditional_theorem",
            readiness="conditional_on_open_input",
            claim="All-shift lambda=-100 entry would complete contiguous order-five forward invariance through zero.",
            formula=conclusions["conditional_forward"],
            proof_boundary="Conditional on the explicit entry target only.",
        ),
        ReductionRow(
            id="co5utfr_11_lower_layer_kill_gate",
            role="non_promotion_guard",
            readiness="ready_to_apply",
            claim="The completed signed layers through order four do not imply order five abstractly.",
            formula="H_(5,0)=-1/3657830400000<0",
            proof_boundary="Exact rational countermodel, not the Xi sequence.",
            diagnostics=countermodel,
        ),
        ReductionRow(
            id="co5utfr_12_open_entry",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove strict all-shift order-five entry at lambda=-100, equivalently strict log-concavity of H_4 there.",
            formula=conclusions["open_entry_target"],
            proof_boundary="This is the sole open input in this order-five propagation route.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_uniform_tail_flow_reduction",
        "date": "2026-07-13",
        "status": (
            "exact uniform eventual order-five tail and conditional forward reduction "
            "with one open lambda=-100 entry"
        ),
        "proof_boundary": (
            "This artifact proves an unconditional compact-uniform eventual "
            "order-five tail and the exact cooperative propagation reduction. "
            "It does not prove all-shift order-five entry, full order-five "
            "invariance, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md",
            "outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md",
            "outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md",
            "outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md",
            "https://arxiv.org/abs/1910.01227",
            "https://arxiv.org/abs/2007.13582",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order5_uniform_tail_flow_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order5_uniform_tail_flow_reduction.py"
        ),
        "source_contract": sources,
        "uniform_suitability_extension": suitability,
        "lambert_extension": lambert,
        "newton_transfer": newton,
        "determinant_cancellation": cancellation,
        "exact_flow": flow,
        "countermodel": countermodel,
        "conclusions": conclusions,
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": sum(
                row.readiness == "ready_to_apply" for row in rows
            ),
            "conditional_rows": sum(
                row.readiness == "conditional_on_open_input" for row in rows
            ),
            "open_rows": sum(
                row.readiness == "not_ready_to_apply" for row in rows
            ),
            "suitability_coefficients_checked": len(suitability["coefficient_checks"]),
            "lambert_derivative_orders_checked": len(lambert["rows"]),
            "newton_coefficients_checked": len(newton["coefficient_rows"]),
            "determinant_permutations_checked": cancellation["permutations_checked"],
            "uniform_eventual_tail_theorems": 1,
            "cooperative_flow_theorems": 1,
            "conditional_forward_theorems": 1,
            "lower_layer_countermodels": 1,
            "open_entry_targets": 1,
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    cancellation = artifact["determinant_cancellation"]
    flow = artifact["exact_flow"]
    conclusions = artifact["conclusions"]
    countermodel = artifact["countermodel"]
    lines = [
        "# Jensen-Window PF Compound Order-Five Uniform Tail And Flow Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact uniform eventual order-five tail and conditional forward",
        "reduction with one open `lambda=-100` entry. This is not a proof of",
        "all-shift order five, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order5_uniform_tail_flow_reduction`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_uniform_tail_flow_reduction.py",
        "```",
        "",
        "## Compact-Heat Asymptotic",
        "",
        "O'Sullivan's all-order suitable-multiplier theorem applies to the same",
        "compact family through every fixed order needed here. The exact",
        "suitability audit and Lambert recurrence now run through order eleven:",
        "",
        "```text",
        conclusions["heat_tilt_extension"],
        conclusions["higher_theta"],
        "```",
        "",
        "Nine-node Newton interpolation transfers those estimates to the graded",
        "ratio coefficients. For `M=n+8` and `h=Delta(M)^2`, the exact",
        "120-permutation determinant expansion is",
        "",
        "```text",
        "[h^0,...,h^10] det K=[0,0,0,0,0,0,0,0,0,0,294912*G_2^10]",
        cancellation["first_nonzero_term"],
        "```",
        "",
        "The coefficient is the positive Vandermonde constant",
        "`2^10*(1!*2!*3!*4!)=294912`; `G_3,...,G_11` cancel from the first",
        "nonzero term. Since `G_2->1` uniformly,",
        "",
        "```text",
        conclusions["uniform_eventual_tail"],
        "```",
        "",
        "## Exact Entry Coordinate",
        "",
        "Desnanot-Jacobi gives",
        "",
        "```text",
        flow["condensation"],
        "```",
        "",
        "The completed signs `H_(3,n)<0` and `H_(4,n)>0` make the order-five",
        "entry target exactly",
        "",
        "```text",
        "H_(4,n+1)(-100)^2>H_(4,n)(-100)*H_(4,n+2)(-100), every n>=0.",
        "```",
        "",
        "Thus the new analytic problem is strict log-concavity of the actual",
        "positive order-four determinant sequence at `lambda=-100`.",
        "",
        "## Cooperative Flow",
        "",
        "Exact affine-Hankel differentiation and the adjacent Plucker identity",
        "give",
        "",
        "```text",
        flow["cooperative_flow"],
        "```",
        "",
        "The off-diagonal coefficient is positive because order four is already",
        "complete. The flow needs no order-six sign. If the open entry target is",
        "proved, the uniform tail confines the remaining system to finitely many",
        "indices, and variation of constants gives",
        "",
        "```text",
        conclusions["conditional_forward"],
        "```",
        "",
        "## Lower-Layer Kill Gate",
        "",
        "Order five is not a formal consequence of the lower signed layers. The",
        "positive rational sequence",
        "",
        "```text",
        str(countermodel["sequence"]),
        "```",
        "",
        "has every available signed contiguous minor through order four strict,",
        "but",
        "",
        "```text",
        f"H_(5,0)={countermodel['H5_n0']}<0,",
        f"H_(4,1)^2-H_(4,0)*H_(4,2)={countermodel['H4_log_concavity_margin']}<0.",
        "```",
        "",
        "## Open Input",
        "",
        "```text",
        "prove H_(5,n)(-100)>0 for every integer n>=0",
        "equivalently prove strict log-concavity of H_(4,n)(-100) for every n>=0",
        "```",
        "",
        "This is now the sole order-five propagation input. It is not supplied by",
        "finite evidence, the lower layers, or the eventual tail.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md",
        "outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
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
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-five uniform tail and flow reduction: "
        f"{summary['rows']} rows, "
        f"{summary['determinant_permutations_checked']} determinant permutations, "
        f"{summary['open_entry_targets']} open entry target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
